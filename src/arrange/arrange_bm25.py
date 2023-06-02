from typing import List

from src.utils.OSM import OSM
from src.utils.Oligo import Oligo


def check_result(best_matches: List[OSM], token="UCGA"):
    if len(best_matches) == 0:
        return 0

    return len([b for b in best_matches if
                (b.oligo.sequence in b.spectrum['params']['title'] or (
                        b.oligo.sequence[::-1] in b.spectrum['params']['title'])) and b.matched_score > 0])
    # return len([b for b in best_matches if b.oligo.sequence in token and b.matched_score > 0])


def check_precursor_filter(oligos: List[Oligo], token="UCGA"):
    if len(oligos) == 0:
        return 0

    return len([o for o in oligos if
                (o.sequence in token or (
                        o.sequence[::-1] in token))])
    # return len([b for b in best_matches if b.oligo.sequence in token and b.matched_score > 0])
