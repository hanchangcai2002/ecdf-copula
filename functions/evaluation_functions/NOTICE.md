# Third-party code notice

The `EVALUATION FUNCTIONS/` and `PREPROCESSING/` folders in this directory are
vendored, largely unmodified, from:

> Hernandez, M., Epelde, G., Alberdi, A., Cilla, R., Rankin, D. (2023).
> Synthetic Tabular Data Evaluation in the Health Domain Covering Resemblance,
> Utility, and Privacy Dimensions. *Methods of Information in Medicine*, 62,
> e19-e38. https://doi.org/10.1055/s-0042-1760247
>
> Source repository: https://github.com/Vicomtech/STDG-evaluation-metrics

Licensed under the MIT License (see `LICENSE.txt`, copyright Open Health
Imaging Foundation). We reuse these evaluation functions (univariate/
multivariate/dimensional resemblance, data labelling, utility, similarity,
membership inference, attribute inference) to evaluate our own eCDF-copula
synthetic data against the same resemblance/utility/privacy framework, as
described in Section 3.3 ("Evaluation Framework") of our paper.

This vendored copy excludes the original repository's own example datasets,
notebooks, and generated results — only the reusable Python function modules
are kept.
