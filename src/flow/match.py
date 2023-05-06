import concurrent.futures
import csv
from typing import List

from src.utils.OSM import OSM
from src.utils.Oligo import Oligo


def cal_delta_ppm(x, y):
    return 1e6 * abs(x - y) / x


def filter_by_oligo_mass(oligos, mass, delta_threshold):
    return [o for o in oligos if abs(cal_delta_ppm(mass[0], o.mass)) < delta_threshold]


def norm_intensities(intensities):
    max_in = max(intensities)
    return intensities / max_in


def cal_bm25(intensities, oligo, peak_matches_info, avg_og_len, K, t):
    # TODO: 参数也许可以优化
    score = 0.0
    len_oligo = len(oligo.sequence)
    norm_ins = norm_intensities(intensities)

    for i, j, delta in peak_matches_info:
        score += (norm_ins[i] * (K + 1) / (
                norm_ins[i] + (K * (1 - t + t * (len_oligo / avg_og_len)))) / delta)

    return score


def cal_norm_cover(intensities, peak_matches_info):
    return len(peak_matches_info) / len(intensities)


def cal_norm_cover_intensities(intensities, peak_matches_info):
    score = 0.0
    max_in = max(intensities)
    for i, j, delta in peak_matches_info:
        score += intensities[i] / max_in
    return score / len(intensities)


def cal_osm(spectrum, oli, args_dict, avg_length):
    sorted_mz_array = spectrum['m/z array']
    peak_matches_info = []
    i = 0
    j = 0
    while i < len(sorted_mz_array) and j < len(oli.masses):
        delta_ppm = cal_delta_ppm(sorted_mz_array[i], oli.masses[j][1])
        if delta_ppm < float(args_dict['peak_window_ppm']):
            peak_matches_info.append((i, j, delta_ppm))
            i += 1
            j += 1
        elif sorted_mz_array[i] < oli.masses[j][1]:
            i += 1
        else:
            j += 1
    # for i in range(len(sorted_mz_array)):
    #     for j in range(len(oli.masses)):
    #         delta_ppm = cal_delta_ppm(sorted_mz_array[i], oli.masses[j][1])
    #         if delta_ppm < float(args_dict['peak_window_ppm']):
    #             peak_matches_info.append((i, j, delta_ppm))
    score = 0
    if avg_length != 0 and len(peak_matches_info) != 0:
        # # norm_intensity
        # score = cal_norm_cover_intensities(spectrum['intensity array'], peak_matches_info)

        # # len(peak_matches_info) / len(intensities)
        # score = cal_norm_cover(spectrum['intensity array'], peak_matches_info)

        # BM25
        score = cal_bm25(spectrum['intensity array'], oli, peak_matches_info, avg_length,
                         float(args_dict['bm25_k']),
                         float(args_dict['bm25_b']))

    with open("peak_match_infos_" + spectrum['params']['title'] + '_' + oli.sequence+'.csv', 'w', newline='') as file:
        file.write("spec_mz_index,oli_mass_index,delta_ppm\n")
        for ps in peak_matches_info:
            file.write("{},{},{}\n".format(ps[0],ps[1],ps[2]))

    return OSM(score, oli, spectrum, peak_matches_info)


def match(spectrum, oligos: List[Oligo], args_dict):
    """

    :param spectrum:
    :param oligos:
    :param args_dict: cfg_dict['match']
    :return:
    """
    # sorted_mz_array = spectrum['m/z array']

    filtered_oligos = filter_by_oligo_mass(oligos, spectrum['params']['pepmass'], int(args_dict['precursor_delta_ppm']))
    if len(filtered_oligos) == 0:
        return OSM(-1, oligos[0], spectrum, [])
    else:
        avg_length = sum([len(o.sequence) for o in filtered_oligos]) / len(filtered_oligos)
        filtered_oligos.sort(key=lambda o: o.mass)

        matches_infos = []
        osm_futures = []

        with concurrent.futures.ThreadPoolExecutor() as executor:
            for oli in filtered_oligos:
                osm_futures.append(executor.submit(cal_osm, spectrum, oli, args_dict, avg_length))

            for future in osm_futures:
                mf = future.result()
                if mf:
                    matches_infos.append(mf)
            # for oli in filtered_oligos:
            #     peak_matches_info = []
            #     i = 0
            #     j = 0
            #     while i < len(sorted_mz_array) and j < len(oli.masses):
            #         delta_ppm = cal_delta_ppm(sorted_mz_array[i], oli.masses[j][1])
            #         if delta_ppm < float(args_dict['peak_window_ppm']):
            #             peak_matches_info.append((i, j, delta_ppm))
            #             i += 1
            #             j += 1
            #         elif sorted_mz_array[i] < oli.masses[j][1]:
            #             i += 1
            #         else:
            #             j += 1
            #     # for i in range(len(sorted_mz_array)):
            #     #     for j in range(len(oli.masses)):
            #     #         delta_ppm = cal_delta_ppm(sorted_mz_array[i], oli.masses[j][1])
            #     #         if delta_ppm < float(args_dict['peak_window_ppm']):
            #     #             peak_matches_info.append((i, j, delta_ppm))
            #     score = 0
            #     if avg_length != 0 and len(peak_matches_info) != 0:
            #         # # norm_intensity
            #         # score = cal_norm_cover_intensities(spectrum['intensity array'], peak_matches_info)
            #
            #         # # len(peak_matches_info) / len(intensities)
            #         # score = cal_norm_cover(spectrum['intensity array'], peak_matches_info)
            #
            #         # BM25
            #         score = cal_bm25(spectrum['intensity array'], oli, peak_matches_info, avg_length,
            #                          float(args_dict['bm25_k']),
            #                          float(args_dict['bm25_b']))
            #     matches_infos.append(OSM(score, oli, spectrum, peak_matches_info))



        best_match = max(matches_infos, key=lambda x: x.matched_score)
        print(best_match.to_string())

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


if __name__ == '__main__':
    print(cal_delta_ppm(100, 100.005))
