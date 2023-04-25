from typing import List

from example.utils.OSM import OSM
from example.utils.Oligo import Oligo


def filter_by_oligo_mass(oligos, mass, delta):
    return [o for o in oligos if abs(o.mass - mass[0]) < delta]


def cal_bm25(intensities, oligo, peak_matches_info, avg_og_len, K, t):
    # TODO: 参数也许可以优化
    score = 0.0
    len_oligo = len(oligo.sequence)

    for i, j, delta in peak_matches_info:
        score += (intensities[i] * (K + 1) / (
                intensities[i] + (K * (1 - t + t * (len_oligo / avg_og_len)))) / delta)

    return score


def match(spectrum, oligos: List[Oligo], args_dict):
    """

    :param spectrum:
    :param oligos:
    :param args_dict: cfg_dict['match']
    :return:
    """
    mz_array = spectrum['m/z array']

    filtered_oligos = filter_by_oligo_mass(oligos, spectrum['params']['pepmass'], int(args_dict['precursor_delta_ppm']))
    avg_length = sum([len(o.sequence) for o in filtered_oligos])

    matches_infos = []
    for oli in filtered_oligos:
        peak_matches_info = []
        for i in range(len(mz_array)):
            for j in range(len(oli.masses)):
                delta = (oli.masses[j][1] - mz_array[i]) / mz_array[i]
                if abs(delta) < int(args_dict['peak_window_ppm']):
                    peak_matches_info.append((i, j, delta))
        score = cal_bm25(spectrum['intensity array'], oli, peak_matches_info, avg_length, float(args_dict['bm25_k']),
                         float(args_dict['bm25_t']))
        matches_infos.append(OSM(score, oli, spectrum, peak_matches_info))

    best_match = max(matches_infos, key=lambda x: x[0])

    return best_match


def get_rna_matches(best_matches: List[OSM]):
    olis = [b.oligo for b in best_matches]
    protein_match = {}
    for oli in olis:
        for loca in oli.sequence_location:
            if loca['molecule'] in protein_match:
                protein_match[loca['molecule']].append(loca)
            else:
                protein_match[loca['molecule']] = [loca]

    return protein_match
