import concurrent.futures
import glob
import os

import pandas as pd
from pyteomics import fasta
from typing import List

from src.RNA_PTM import read_fasta, prepare_for_silicon_spec
from src.flow.mod_digest import mod_0_1_2_mode
from src.utils.Break import Break
from src.utils.Enzyme import create_enzymes_from_table, get_used_enzymes
from src.utils.Oligo import Oligo
from src.utils.cfg import read_cfg, get_cfg_items


def find_cut_poses(rna, end5, end3, used_enzymes):
    # sites = []
    cut_poses = [{'cut_index': 0, 'cut_end5': "-1", 'cut_end3': end5, 'enzyme_name': "seq_start"}]

    seq = rna[1]
    for ezy in used_enzymes:
        index = -1
        for cut_info in used_enzymes[ezy].cut_list:
            while True:
                index = seq.find(cut_info[0], index + 1)
                if index == -1:
                    break
                cut_poses.append({'cut_index': index + cut_info[1], 'cut_end5': used_enzymes[ezy].cut_5end,
                                  'cut_end3': used_enzymes[ezy].cut_3end, 'enzyme_name': ezy})
    cut_poses.append({'cut_index': len(seq), 'cut_end5': end3, 'cut_end3': "-1", 'enzyme_name': "seq_end"})
    cut_poses = sorted(cut_poses, key=lambda x: x['cut_index'])
    return cut_poses


def without_digest(rna, enzymes_h, args, is_decoy=False):
    return [Oligo(
        sequence=rna[1],
        miss=[0],
        chemistry_5=args["random_cut_end3"],
        chemistry_3=args["random_cut_end5"],
        sequence_location=[{'molecule': rna[0],
                            'residue_start': 0,
                            'residue_end': len(rna[1]),
                            'enzyme_name_start': "without_digest",
                            'enzyme_name_end': "without_digest",
                            'miss': 0,
                            'mod_infos': [],
                            'is_decoy': is_decoy}]
    )]


def specific(rna, enzymes_s, args, is_decoy=False):
    seq: str = rna[1]
    min_len = int(args["min_oligo_length"])
    max_len = int(args["max_oligo_length"])
    miss = int(args["max_miss_site"])

    cut_poses = find_cut_poses(rna, args['end5'], args['end3'], enzymes_s)
    # for cut_pos in cut_poses:
    oligos = []
    try:
        for start_p in range(len(cut_poses) - 1):
            for end_p in range(start_p + 1, len(cut_poses)):
                if (end_p - start_p) > (miss + 1):
                    raise Break()
                oligo_length = cut_poses[end_p]['cut_index'] - cut_poses[start_p]['cut_index']
                if max_len >= oligo_length >= min_len:
                    seq_to_add = seq[cut_poses[start_p]['cut_index']:cut_poses[end_p]['cut_index']]
                    # print(seq_to_add))
                    oligo = Oligo(
                        sequence=seq_to_add,
                        miss=[miss],
                        chemistry_5=cut_poses[start_p]['cut_end3'],
                        chemistry_3=cut_poses[end_p]['cut_end5'],
                        sequence_location=[{'molecule': rna[0], 'residue_start': cut_poses[start_p]['cut_index'],
                                            'residue_end': cut_poses[end_p]['cut_index'],
                                            'enzyme_name_start': cut_poses[start_p]['enzyme_name'],
                                            'enzyme_name_end': cut_poses[end_p]['enzyme_name'],
                                            'mod_infos': [],
                                            'is_decoy': is_decoy}]
                    )
                    oligos.append(oligo)
                    # print(oligo.to_string())
    except Break:
        pass

    return oligos


def half_specific(rna, enzymes_h, args, is_decoy=False):
    seq: str = rna[1]
    min_len = int(args["min_oligo_length"])
    max_len = int(args["max_oligo_length"])
    miss = int(args["max_miss_site"])

    cut_poses = find_cut_poses(rna, args['end5'], args['end3'], enzymes_h)

    oligos = []
    # 从前往后遍历，即前特后非
    for start_p in range(1, len(cut_poses) - 1):
        for i_len in range(min_len, max_len + 1):
            end_pos = cut_poses[start_p]['cut_index'] + i_len
            if end_pos > cut_poses[start_p + miss + 1]['cut_index']:
                break
            if end_pos <= len(seq):
                seq_to_add = seq[cut_poses[start_p]['cut_index']: end_pos]
                # print(seq_to_add))
                oligo = Oligo(
                    sequence=seq_to_add,
                    miss=[miss],
                    chemistry_5=cut_poses[start_p]['cut_end3'],
                    chemistry_3=args["random_cut_end5"],
                    sequence_location=[{'molecule': rna[0], 'residue_start': cut_poses[start_p]['cut_index'],
                                        'residue_end': end_pos,
                                        'enzyme_name_start': cut_poses[start_p]['enzyme_name'],
                                        'enzyme_name_end': "half_random",
                                        'miss': miss,
                                        'mod_infos': [],
                                        'is_decoy': is_decoy}]
                )
                oligos.append(oligo)
                # print(oligo.to_string())

    # 从后往前遍历，即前非后特
    for end_p in range(len(cut_poses) - 2, -1, -1):
        for i_len in range(min_len, max_len + 1):
            start_pos = cut_poses[end_p]['cut_index'] - i_len
            if cut_poses[end_p - miss - 1]['cut_index'] > start_pos:
                break
            if start_pos >= 0:
                seq_to_add = seq[start_pos: cut_poses[end_p]['cut_index']]
                # print(seq_to_add))
                oligo = Oligo(
                    sequence=seq_to_add,
                    miss=0,
                    chemistry_5=args["random_cut_end3"],
                    chemistry_3=cut_poses[end_p]['cut_end5'],
                    sequence_location=[{'molecule': rna[0], 'residue_start': start_pos,
                                        'residue_end': cut_poses[end_p]['cut_index'],
                                        'enzyme_name_start': "half_random",
                                        'enzyme_name_end': cut_poses[end_p]['enzyme_name'],
                                        'miss': 0,
                                        'mod_infos': [],
                                        'is_decoy': is_decoy}]
                )
                oligos.append(oligo)
                # print(oligo.to_string())
    return oligos


def nonspecific(rna, enzymes_h, args, is_decoy=False):
    seq: str = rna[1]
    min_len = int(args["min_oligo_length"])
    max_len = int(args["max_oligo_length"])

    oligos = []
    for i in range(min_len, max_len):
        if i <= len(seq):
            for pos in range(len(seq) - i):
                seq_to_add = seq[pos: pos + i + 1]
                # print(seq_to_add))
                oligo = Oligo(
                    sequence=seq_to_add,
                    miss=[0],
                    chemistry_5=args["random_cut_end3"],
                    chemistry_3=args["random_cut_end5"],
                    sequence_location=[{'molecule': rna[0],
                                        'residue_start': pos,
                                        'residue_end': pos + i,
                                        'enzyme_name_start': "non_random",
                                        'enzyme_name_end': "non_random",
                                        'miss': 0,
                                        'mod_infos': [],
                                        'is_decoy': is_decoy}]
                )
                oligos.append(oligo)
                # print(oligo.to_string())
    return oligos


def shrink_oligos(oligos: List[Oligo]):
    """
    基于sequence、chemistry_5、chemistry_3、修饰位置对Oligos列表去冗余
    只有序列、修饰都完全一样的oligo才会被合并
    :param oligos:
    :return:
    """
    # LOG: 对oligos按照序列与左右end去冗余（对已有的list与冗余效率不知道是否足够高--》方便上一步进行多线程）
    new_list: List[Oligo] = []
    seen_values: list = []

    for i in range(len(oligos)):
        oligo = oligos[i][0]
        print(oligo)

        value = (oligo.sequence, oligo.chemistry_5, oligo.chemistry_3, oligo.sequence_location[0]['mod_infos'])
        index = seen_values.index(value)
        if index != -1:
            new_list[index].sequence_location.append(oligo.sequence_location)
        else:
            new_list.append(oligo)
            seen_values.append(value)

    print(new_list)
    return new_list


def digest_to_oligos(args, rna, enzymes, is_decoy=False):
    enzyme_type = int(args["enzyme_type"])
    olgs = []
    # for rna in rnas:
    if enzyme_type == -1:
        olgs = without_digest(rna, enzymes, args, is_decoy)
    elif enzyme_type == 0:
        olgs = specific(rna, enzymes, args, is_decoy)
    elif enzyme_type == 1:
        olgs = half_specific(rna, enzymes, args, is_decoy)
    elif enzyme_type == 2:
        olgs = nonspecific(rna, enzymes, args, is_decoy)

    return olgs

# if __name__ == '__main__':
#     # 读取cfg文件
#     config = read_cfg('data/cfg.ini')
#     paths_dict = get_cfg_items(config, "paths")
#
#     # FASTA文件处理
#     fasta_datas = read_fasta(paths_dict['fasta_dir'])
#
#     args_dict = get_cfg_items(config, "digest")
#
#     used_enzymes = get_used_enzymes(paths_dict['enzyme_path'], args_dict['enzymes_name'])
#
#     all_oligos = []
#     futures = []
#     with concurrent.futures.ThreadPoolExecutor() as executor:
#         for f in fasta_datas:
#             while True:
#                 try:
#                     header, sequence = next(f)
#                     #
#                     # oligos = digest_to_oligos(args_dict, (header, sequence), used_enzymes)
#                     futures.append(executor.submit(digest_to_oligos, args_dict, (header, sequence), used_enzymes))
#                     # Print the header and sequence of the first protein
#                     print('Header:', header)
#                     print('Sequence:', sequence)
#                 except StopIteration:
#                     # 当没有更多的RNA序列时，停止循环
#                     break
#
#         # 将所有的oligo整合，然后去冗余
#         # TODO 所有oligo在集成在一个表里面，似乎会占用大量空间，后续应考虑是否需优化
#         for future in futures:
#             all_oligos.append(future.result())
#
#     # all_oligos = shrink_oligos(all_oligos)
#     # # all_oligos = shrink_oligos(all_oligos)
#     elements, nucleotides, mods = prepare_for_silicon_spec()
#     mod_0_1_2_mode(paths_dict['known_mod_path'], mods, all_oligos)
