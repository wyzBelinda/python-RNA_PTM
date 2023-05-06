from typing import List

from src.utils.OSM import OSM


def check_result(best_matches: List[OSM], token="AGUC"):
    if len(best_matches) == 0:
        return 0
    token = best_matches[0].oligo.sequence_location[0]['molecule']
    return len([b for b in best_matches if
                b.oligo.sequence in b.spectrum['params']['title'] and b.matched_score > 0])
    # return len([b for b in best_matches if b.oligo.sequence == token and b.matched_score > 0])
