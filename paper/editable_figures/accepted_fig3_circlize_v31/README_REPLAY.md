# Figure 3a replay

This is the actual R/circlize/ggseqlogo adaptation of the supplied plotting code. Original code and reference PDF are preserved under `provided_reference_code`; CARBON inputs remain separate under `inputs`.

On Windows, copy this directory to an ASCII-only working directory. The verified runtime was `D:/Anaconda/envs/repro_sctype/Scripts/Rscript.exe` (R 4.3.3). Add the environment root, `Library/bin` and `Scripts` to PATH; set LANG and LC_ALL to C. From that working directory run:

```
Rscript install_offline.R
Rscript render_circos.R 285 circos_carbon_285mm
```

The isolated `Rlibrary` is used; CRAN source packages and their hashes are bundled. Existing environment packages such as ggplot2, gridExtra and colorspace are described in `R_SESSION_RENDER.txt`. The original online installer at revision root was an unsuccessful first attempt; this bundled offline installer is the replay entry.

The renderer writes a vector PDF and geometry CSV. The final intact 85 x 98.5 mm Figure 3 and five-page S21 are also supplied. `compose_and_package.py` records their construction in the original revision directory and clips b/c from the preserved v30 PDF without modifying those panels. The current editable master is R plus native PDF; old v30 SVG/PNG/panel-a aliases are archived separately in `manuscript/figures/previous_v30_fig3_formats`.

All 219 tasks, 150 RBP symbols, frozen 9-nt PPM crops, Euclidean distances, average-linkage tree, heights and leaf order remain unchanged. Six pastel cutree(k=6) sectors are visual orientation only (196/11/4/6/1/1), not biological categories. No paired/unpaired RNA structure channels are inferred. Full names are tiny in the single-column overview and provided in the enlarged export and S21. R emitted ggseqlogo deprecation/rotated-viewport warnings; rendered outputs were visually checked. No model inference or significance tests were rerun.
