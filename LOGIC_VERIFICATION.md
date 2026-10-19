# Checked demonstration and limits

Run python scripts/verify_demo_cpu.py from the repository root. The shell compatibility entry scripts/verify_logic_regression.sh calls the same verifier. It loads the real included AARS checkpoint, explicitly selects CPU with --CUDA -1, passes --task, evaluates 1,283 held-out demo windows and compares six metrics against the saved demo result with absolute tolerance 1e-5. Reports go to outputs/demo_cpu.

This checks the included demonstration. It does not prove equality between the current CARBON implementation and historical MARCP code, fresh training reproducibility, fair three-model scoring, or all-paper reproduction. The broken historical AGGF1 regression command has been replaced; no equivalence result is fabricated.

The demo chr1/2/3 split is an example only. Formal archived experiment manifests use train chromosomes 3,4,5,6,7,10,11,12,13,14,17,18,19,20,21,22; validation 2,9,16; test 1,8,15. Follow each experiment's own saved configuration.
