import concurrent.futures
import glob
import os.path

import pandas as pd
from pyteomics import mgf, fasta

from example.package_2.filter_FDR import tda_fdr_test
from example.package_2.match import match, get_rna_matches
from example.package_2.mod_digest import find_mods
from example.utils.Element import create_elements_from_file
from example.utils.Enzyme import get_used_enzymes
from example.utils.Modification import create_modifications_from_table, get_known_mod
from example.utils.Nucleotide import create_nucleotides_from_table
from example.utils.cfg import read_cfg, get_cfg_items


def prepare_for_silicon_spec():
    """
    运行
    :return:
    """
    # 创建 element 对象列表
    elements = create_elements_from_file('data/element.ini')
    # elements = create_elements_from_file(paths_dict['element_path'])

    # 读入表格
    df = pd.read_csv('data/nts.ini')

    # 创建 nucleotide 对象列表
    nucleotides = create_nucleotides_from_table(df)
    # 计算每个基团的质量并赋值
    for nt in nucleotides:
        nucleotides[nt].gen_mass_by_elements(elements)

    # 载入所有已知的修饰数据
    # 读入表格
    df = pd.read_csv('data/mods.csv')

    # 创建 nucleotide 对象列表
    mods = create_modifications_from_table(df)

    return elements, nucleotides, mods


def read_fasta(fasta_path):
    # FASTA文件处理
    if os.path.isdir(fasta_path):
        # 获取指定目录中的所有FASTA文件路径
        fasta_paths = glob.glob(os.path.join(fasta_path, '*.fasta'))
    else:
        fasta_paths = [fasta_path]

    fasta_datas = []
    for path in fasta_paths:
        fasta_datas.append(fasta.read(path))
    # with fasta.read('data/for_test_small.fasta') as f:

    return fasta_datas


def read_mgf(mgf_path):
    # FASTA文件处理
    if os.path.isdir(mgf_path):
        # 获取指定目录中的所有FASTA文件路径
        mgf_paths = glob.glob(os.path.join(mgf_path, '*.mgf'))
    else:
        mgf_paths = [mgf_path]

    mgf_datas = []
    for path in mgf_paths:
        mgf_datas.append(mgf.read(path))
    # with fasta.read('data/for_test_small.fasta') as f:

    return mgf_datas


def gen_rev2_decoy(rna):
    return "REV2_" + rna[0], rna[1]


if __name__ == '__main__':

    # 读取cfg文件
    config = read_cfg('data/cfg.ini')
    paths_dict = get_cfg_items(config, "paths")

    # FASTA文件处理
    mgf_datas = read_fasta(paths_dict['fasta_dir'])

    # -----酶切------
    digest_dict = get_cfg_items(config, "digest")

    used_enzymes = get_used_enzymes(paths_dict['enzyme_path'], digest_dict['enzymes_name'])

    all_oligos = []
    futures = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for f in mgf_datas:
            while True:
                try:
                    header, sequence = next(f)
                    from example.package_2.enzyme_digest import digest_to_oligos, shrink_oligos

                    # oligos = digest_to_oligos(args_dict, (header, sequence), used_enzymes)
                    futures.append(executor.submit(digest_to_oligos, digest_dict, (header, sequence), used_enzymes))
                    futures.append(
                        executor.submit(digest_to_oligos, digest_dict, (gen_rev2_decoy((header, sequence))),
                                        used_enzymes, True))

                    # Print the header and sequence of the first protein
                    print('Header:', header)
                    print('Sequence:', sequence)
                except StopIteration:
                    # 当没有更多的RNA序列时，停止循环
                    break

        # 将所有的oligo整合，然后去冗余
        # TODO 所有oligo在集成在一个表里面，似乎会占用大量空间，后续应考虑是否需优化
        for future in futures:
            o = future.result()
            if o:
                all_oligos.extend(o)

    match_dict = get_cfg_items(config, "match")

    elements, nucleotides, mods = prepare_for_silicon_spec()

    known_mod = get_known_mod(paths_dict['known_mod_path'])
    all_oligos = find_mods(known_mod, all_oligos)
    # all_oligos=mod_0_1_2_mode(paths_dict['known_mod_path'], mods, all_oligos)

    # all_oligos = shrink_oligos(all_oligos)
    for oli in all_oligos:
        oli.cal_fragments_mass(elements, nucleotides, mods)

    mgf_datas = read_mgf(paths_dict['mgf_dir'])
    best_matches = []
    match_futures = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for m in mgf_datas:
            for spectrum in m:
                print(spectrum['params'])
                # print(spectrum['m/z array'])
                # print(spectrum['intensity array'])
                futures.append(executor.submit(match, spectrum, all_oligos, match_dict))

        for future in match_futures:
            mf = future.result()
            if mf:
                best_matches.extend(mf)

    protein_match = get_rna_matches(best_matches)
    print(protein_match)

    # -------------
    # 分析FDR&过滤
    analysis_dict = get_cfg_items(config, "analysis")
    n_target, n_decoy, n_pass = tda_fdr_test(best_matches, analysis_dict['fdr_threshold'])

    # -------------
    # 可视化

    # # ---------------------
    # mgf_dir = "data/"
    # # 获取指定目录中的所有mgf文件路径
    # mgf_paths = glob.glob(os.path.join(mgf_dir, '*.mgf'))
    #
    # for mgf_path in mgf_paths:
    #     # Open the MGF file and read the whole spectrum
    #     with mgf.read(mgf_path) as spectra:
    #         # with mgf.read('data/02_OH-AGUC-OH_CID_srx.mgf') as spectra:
    #         while True:
    #             try:
    #                 spectrum = next(spectra)
    #                 print(spectrum['TITLE'])
    #
    #                 # # Print the m/z and intensity values for each peak in the first spectrum
    #                 # for mz, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):
    #                 #     print('m/z:', mz, 'Intensity:', intensity)
    #             except StopIteration:
    #                 # 当没有更多的谱图时，停止循环
    #                 break
