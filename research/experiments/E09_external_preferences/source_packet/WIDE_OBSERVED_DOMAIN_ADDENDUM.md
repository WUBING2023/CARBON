# Missing native word associations — scope repair before successful group evaluation

The complete controlled-response and external-label axes contain all1,024 five-mers. However,1,024 native context windows for some proteins do not contain every possible word. A missing native word cannot be assigned a fabricated Phi association of zero.

The primary full-word evaluation retains all1,024 words for Full,Sequence,frozenFlat,frozenPerm,native frequency and all four universal base-count predictors. Native frequency zero is a real absence count. Native Phi–word association is moved to an explicitly auxiliary shared observed-word domain. For each task, all model and simple predictors are recomputed on exactly the same words with at least one native occurrence; same-composition pairs also require both words in that domain. Every task's observed/missing word count is exported. No Full full-domain rho is compared to Phi partial-domain rho.

In the controlled substitution design the genuine Phi-only input is identical across all words, because the original conservation vector is frozen. Its word-ordering prediction is therefore tied: rho is undefined/NA and same-composition concordance is0.5. This is an input-only logical baseline, not a newly fitted prediction model. Its NA rho is never replaced by zero or used in an apparent mean difference.

This correction was prompted by a source-availability assertion before the first successful complete28-protein evaluation. It does not change a model checkpoint, response, k-mer axis, source label, experiment inclusion, primary/secondary endpoint or inference context. All extra native-Phi comparisons are labeled with their actual common domain.
