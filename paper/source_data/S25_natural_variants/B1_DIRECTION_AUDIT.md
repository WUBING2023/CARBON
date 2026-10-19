# Verified QKI contrast direction

The byte-exact accepted B1 input packet and its preparation/render scripts are copied into `b1_reference/`; originals outside this work directory were not edited.

`QKI_words_RNA` is ordered UCAAC, CAUAC, CUACA, CUAAC. Corresponding `QKI_assayZ` values are 0.50133013, 0.45951894, 0.55066383, 2.61200781. The preparation code computes `z[:,-1,None] - z[:,:3]`: CUAAC minus each neighbor. Therefore the positive contrasts and the rendered arrows both represent **UCAAC/CAUAC/CUACA → CUAAC**.

The three assay changes are +2.11068, +2.15249, +2.06134 Z units; Full mean predicted log-count changes are +0.294449, +0.267732, +0.289454. Existing accepted figure arrows are correct. Any prose saying CUAAC → neighbors for these positive numbers must be corrected. The original a/b/c lettering and 85 × 132 mm geometry remain available unchanged.
