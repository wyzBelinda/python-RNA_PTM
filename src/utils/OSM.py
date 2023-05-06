from src.utils.Oligo import Oligo


class OSM:
    def __init__(self, matched_score, oligo: Oligo, spectrum, peak_matches_info, tda_fdr=0):
        self.matched_score = matched_score
        self.oligo = oligo
        self.spectrum = spectrum
        self.peak_matches_info = peak_matches_info
        self.tda_fdr = tda_fdr

    def to_string(self) -> str:
        return f"{self.matched_score} {self.oligo.sequence} {self.oligo.mass} {self.spectrum['params']['pepmass'][0]} {self.tda_fdr} {self.oligo.sequence_location[0]['is_decoy']}"
