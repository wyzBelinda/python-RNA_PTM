import concurrent.futures
import glob
import os.path

import pandas as pd
from pyteomics import mgf, fasta

from src.arrange.arrange_bm25 import check_result
from src.flow.filter_FDR import tda_fdr_test
from src.flow.match import match, get_rna_matches
from src.flow.mod_digest import find_mods, mod_0_1_2_mode
from src.utils.Element import create_elements_from_file
from src.utils.Enzyme import get_used_enzymes
from src.utils.Modification import create_modifications_from_table, get_known_mod, cal_mods_mass
from src.utils.Nucleotide import create_nucleotides_from_table
from src.utils.cfg import read_cfg, get_cfg_items, backup_cfgs


def prepare_for_silicon_spec(args_dict):
    """
    运行
    :return:
    """
    # 创建 element 对象列表
    elements = create_elements_from_file(args_dict['element_path'])
    # elements = create_elements_from_file(paths_dict['element_path'])

    # 读入表格
    df = pd.read_csv(args_dict['icon_path'])

    # 创建 nucleotide 对象列表
    nucleotides = create_nucleotides_from_table(df)
    # 计算每个基团的质量并赋值
    for nt in nucleotides:
        nucleotides[nt].gen_mass_by_elements(elements)

    # 载入所有已知的修饰数据
    # 读入表格
    df = pd.read_csv(args_dict['mods_path'])

    # 创建 nucleotide 对象列表
    mods = create_modifications_from_table(df)
    cal_mods_mass(mods, elements)

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
    return "REV2_" + rna[0], rna[1][:-1][::-1] + rna[1][-1]


# TODO BM25参数测试结束后，将下行删除
def flows_continuing(test_k, test_b):
    # def flows_continuing():
    # 读取cfg文件
    config = read_cfg('data/cfg.ini')
    paths_dict = get_cfg_items(config, "paths")
    if not os.path.exists(paths_dict['out_dir']):
        os.makedirs(paths_dict['out_dir'])

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

                    from src.flow.enzyme_digest import digest_to_oligos
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
    # TODO:BM25参数测试结束后，将下行删除
    match_dict['bm25_k'] = test_k
    match_dict['bm25_b'] = test_b

    elements, nucleotides, mods = prepare_for_silicon_spec(paths_dict)

    known_mod = get_known_mod(paths_dict['known_mod_path'])
    all_oligos = find_mods(known_mod, all_oligos)

    from src.flow.enzyme_digest import shrink_oligos
    all_oligos = shrink_oligos(all_oligos)

    # 碎片离子质量计算
    for oli in all_oligos:
        oli.cal_fragments_mass(elements, nucleotides, mods, int(digest_dict['max_charge']))

        with open(os.path.join(paths_dict['out_dir'], "oligo_ions_" + oli.chemistry_5 + "_" + oli.sequence
                                                      + "_" + oli.chemistry_3 + '.csv'), 'w', newline='') as file:
            file.write("ion_label,ion_mass\n")
            for best_match in oli.masses:
                file.write("{},{}\n".format(best_match[0], best_match[1]))

    mgf_datas = read_mgf(paths_dict['mgf_dir'])
    best_matches = []
    match_futures = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for m in mgf_datas:
            for spectrum in m:
                # print(spectrum['params'])

                match_futures.append(executor.submit(match, spectrum, all_oligos, match_dict, elements['A'][0][0],
                                                     paths_dict['out_dir']))

                with open(os.path.join(paths_dict['out_dir'],
                                       "spectrum_" + spectrum['params']['title'].replace(" ", "_").replace(",",
                                                                                                           "_").replace(
                                           ".", "_").replace('+', '_').replace("/", "_") + ".csv"), 'w',
                          newline='') as file:
                    file.write("mz,intensity\n")
                    for i in range(len(spectrum["m/z array"])):
                        file.write("{},{}\n".format(spectrum["m/z array"][i], spectrum['intensity array'][i]))

        for future in match_futures:
            mf = future.result()
            if mf:
                best_matches.append(mf)

        with open(os.path.join(paths_dict['out_dir'], "best_matches.csv"), 'w', newline='') as file:
            file.write("matched_score,oligo_sequence,oligo_mass,spectrum_title,spectrum_pepmass,tda_fdr,is_decoy\n")
            for best_match in best_matches:
                file.write(best_match.to_string().replace(" ", ",") + "\n")

        print(check_result(best_matches))
        # return check_result(best_matches)

    # protein_match = get_rna_matches(best_matches)
    # print(protein_match)
    #
    # -------------
    # 分析FDR&过滤
    analysis_dict = get_cfg_items(config, "analysis")
    n_target, n_decoy, n_pass = tda_fdr_test(best_matches, float(analysis_dict['fdr_threshold']))
    #
    # # -------------
    # # 可视化

    # 备份参数文件
    backup_cfgs(config['paths']['out_dir'], config)


if __name__ == '__main__':
    # flows_continuing()

    print(flows_continuing(2.1, 2.7))
