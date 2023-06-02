from typing import List

from src.utils.OSM import OSM


def calcFDR(n_target, n_decoy):
    if n_target + n_decoy == 0:
        return 0
    else:
        return n_decoy / n_target


def tda_fdr_test(osm_results: List[OSM], fdr_threshold):
    n_target = 0
    n_decoy = 0
    osm_results.sort(key=lambda x: x.matched_score, reverse=True)
    for result in osm_results:
        if len([d for d in result.oligo.sequence_location if d['is_decoy'] is False]) != 0:
            n_target += 1
        else:
            n_decoy += 1
        result.tda_fdr = calcFDR(n_target, n_decoy)
    osm_results.sort(key=lambda x: x.tda_fdr)
    n_pass = 0
    for result in osm_results:
        if result.tda_fdr <= fdr_threshold:
            result.is_pass = True
            n_pass += 1
        else:
            result.is_pass = False
    print("FDR filter: {} passed, {} filtered.".format(n_pass, len(osm_results) - n_pass))
    return n_target, n_decoy, n_pass
