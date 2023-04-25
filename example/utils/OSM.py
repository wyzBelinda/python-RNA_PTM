class OSM:
    def __init__(self, matched_score, oligo, spectrum, peak_matches_info, tda_fdr=0):
        self.matched_score = matched_score
        self.oligo = oligo
        self.spectrum = spectrum
        self.peak_matches_info = peak_matches_info
        self.tda_fdr = tda_fdr
