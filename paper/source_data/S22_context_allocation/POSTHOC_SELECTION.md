# Posthoc illustrative selection

The case is chosen after Full, NoGRU and near-capacity Local+pos outputs were available. It is not an independent validation set or a preregistered success rate.

Start with all 1,851 original eligible pairs. Require identical local seven-mers, absolute input-adjusted observed share margin at least0.2, at least10 raw eCLIP counts in the two supports, correct Full direction, and wrong directions for both controls. Exclude exact predicted ties. Rank by the smaller of the two control-minus-Full squared-error improvements, descending; break exact-score ties by pair hash ascending. No normalization, exclusion, fit or inferential endpoint changes.

Complete every-row flags, scores and ranks are in SELECTION_ALL_1851.csv. Final candidates are in SELECTION_CANDIDATES.csv; exact funnel, IDs, shares, seeds and native slice checks are in POSTHOC_SELECTION.json. Source shape validation is separate from visual approval.

ZMYM4 is one selected native pair; the atlas gene row averages two eligible pairs. The old LEPR view remains in SUPPLEMENT_PREVIOUS_LEPR.pdf/png/svg with its original source and caption. OCRL failure stays in the main figure.
