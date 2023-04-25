from example.utils.Modification import get_known_mod


def find_mod_poses(oligo, known_mods):
    mods_poses = []
    # [{'mod_index': 0,  'mod_name': "seq_start"}]
    for mod in known_mods:
        known_mods[mod]
        index = -1
        for cut_info in known_mods[mod].cut_list:
            # while True:
            #     index = seq.find(cut_info[0], index + 1)
            #     if index == -1:
            #         break
            mods_poses.append({'cut_index': index + cut_info[1], 'cut_end5': known_mods[mod].cut_5end,
                               'cut_end3': known_mods[mod].cut_3end, 'enzyme_name': mod})
    # mods_poses.append({'cut_index': len(seq), 'cut_end5': end3, 'cut_end3': "-1", 'enzyme_name': "seq_end"})
    cut_poses = sorted(mods_poses, key=lambda x: x['cut_index'])
    return cut_poses


def find_mods(known_mod, all_oligos):
    for known_m in known_mod:
        for i, oligo in enumerate(all_oligos):
            for loca in oligo.sequence_location:
                if loca['residue_start'] <= known_m['Position'] <= loca['residue_end'] \
                        and known_m['Molecule'] == loca['molecule']:
                    # Implement only the modified line if modification is not requested, deleting the
                    # non-modified entry line
                    # if known_m['Include'] == 1:
                    mod_info = {'pos': known_m['Position'] - loca['residue_start'],
                                'mod_symbol': known_m['ID_ext']}
                    s = list(oligo.sequence)
                    if "@" in oligo.mod:

                        s.insert(mod_info['pos'], mod_info['mod_symbol'])
                        # s.insert(known_m['Position'] - loca['residue_start'], known_m['ID_ext'])
                        # s.insert(known_m['Position'] - loca['residue_start'], '[' + known_m['ID_ext'] + ']')

                        oligo.mod = ''.join(s) + "@" * oligo.mod.count('@') + "@"
                        loca['mod_infos'].append(mod_info)
                    else:
                        s.insert(mod_info['pos'], mod_info['mod_symbol'])
                        # s.insert(known_m['Position'] - loca['residue_start'], known_m['ID_ext'])
                        # s.insert(known_m['Position'] - loca['residue_start'], '[' + known_m['ID_ext'] + ']')
                        oligo.mod = oligo.mod + " " + "".join(s) + "@"
                        loca['mod_infos'].append(mod_info)
            print(oligo.to_string())
    return all_oligos

def mod_0_1_2_mode(known_mod_path, mods, all_oligos):
    """
    Add modification lines for nucleotides.
    Copies of the unmodified and modified sequences are appended if requested
    """

    # Add modified bases with the one-letter and extended IDs
    with open(known_mod_path, 'r') as infile:
        for modline in infile:

            if modline.split()[-1].isdigit():
                # Check if the modification is "on" (1 or 2 from the modification file)
                if int(modline.split()[-1]) > 0:
                    mod_lines, positions, del_lines = [], [], []

                    for i, oligo in enumerate(all_oligos):
                        for loca in oligo.sequence_location:
                            if loca['residue_start'] <= int(modline.split()[1]) <= loca['residue_end'] \
                                    and modline.split()[0] == loca['molecule']:
                                # Implement only the modified line if modification is not requested, deleting the
                                # non-modified entry line
                                if modline.split()[-1] == "1":
                                    mod_info = {'pos': int(modline.split()[1]) - loca['residue_start'],
                                                'mod_symbol': modline.split()[2]}

                                    # Tracking trick: add @s at the end of the line for each modified base in the
                                    # fragment
                                    s = list(oligo.sequence)
                                    if "@" in oligo.mod:

                                        s.insert(mod_info['pos'], mod_info['mod_symbol'])

                                        oligo.mod = ''.join(s) + "@" * oligo.mod.count('@') + "@"

                                    else:
                                        s.insert(mod_info['pos'], '[' + mod_info['mod_symbol'] + ']')
                                        oligo.mod = oligo.mod + " " + "".join(s) + "@"

                                    # del_lines.append(line)
                                    # mod_lines.append(" ".join(split)), positions.append(i)

                                # Implement both modified and non-modified lines if requested
                                elif modline.split()[-1] == "2":
                                    mod_info = {'pos': int(modline.split()[1]) - loca['residue_start'],
                                                'mod_symbol': modline.split()[2]}

                                    # Tracking trick: add @s at the end of the line for each modified base in the
                                    # fragment
                                    s = list(oligo.sequence)
                                    if "@" in oligo.mod:
                                        # s = list(oligo.sequence)
                                        s.insert(mod_info['pos'], '[' + mod_info['mod_symbol'] + ']')

                                        oligo.mod = "".join(s) + "@" * oligo.mod.count('@') + "@"

                                    else:
                                        s[int(modline.split()[1]) - loca['residue_start']] = mods[modline.split()[2]]
                                        oligo.mod = oligo.mod + " " + "".join(s) + "@"

                                # Raise an error if an invalid value is specified for the modification option
                                # (value different from 0,1,2)
                                else:
                                    print(
                                        "ERROR! Invalid line {} in the modification file. "
                                        "Execution terminated without output".format(modline))
                                    return  # sys.exit(1)
                        print(oligo.to_string())
    return all_oligos


if __name__ == '__main__':
    known_mod = get_known_mod("data/considered_mods_demo.txt")
    # find_mods(known_mod,all_oligos=[])
