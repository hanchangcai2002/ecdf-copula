# eCDF-Copula: Synthetic Longitudinal Tabular Data Generation via Copula

Code, intermediate data, and evaluation notebooks for:

> Cai H, Yu W, Lu R, Chattopadhyay I, Zhang X, Liu J. *Synthetic Longitudinal
> Tabular Data Generation via Copula.* Biostatistics (in press).

The eCDF-copula method combines the empirical CDF (marginals) with a Gaussian
copula (dependence structure) to synthesize longitudinal tabular health data,
preserving within- and between-visit correlation. It is benchmarked against
Gaussian Multivariate (GM), SDV/PAR, CTGAN, and WGAN-GP on two datasets:

- **Dataset A — REMBRANDT**: n=120, 28 longitudinal variables (MADRS + HARS
  scores across 14 visits). Imputed and synthesized **once**.
- **Dataset B — CHAP**: n=3,612, 20 features (10 binned SBP + 10 binned DBP
  values across pregnancy). Imputed **M=5** times, with **L=10** synthetic
  replicates drawn per imputation (50-dataset ensemble), enabling the
  imputation-vs-synthesis variance decomposition in Figure 7.

## Repository structure

```
functions/            All reusable code: synthesis algorithms, imputation,
                       synthetic-data generation, and evaluation functions.
  synthesis_algorithms/  Core eCDF-copula method (Algorithm 1 & 2 in the paper)
  evaluation_functions/  Vendored resemblance/utility/privacy metric functions
                         (Hernandez et al. 2023, MIT license — see NOTICE.md)
  rembrandt/             Dataset A: MICE imputation + per-method generation
  chap/                  Dataset B: binning + MICE + eCDF-copula/CTGAN pipeline
                         (produces the 50-dataset ensemble) + per-method generation

data/
  raw/                   Real (imputed) data, split by dataset — kept locally
                         only, NOT included in this public repo (see below)
  processed/              Synthetic datasets produced by each method, split by
                         dataset (chap/ also has the 50-run ensemble used for
                         Figure 7)

results/                Evaluation notebooks (.ipynb) and their output
                         metrics/figures — only analyses reported in the paper
  rembrandt/
    resemblance/          Figures 1-2, Table 2 (URA, MRA)
    privacy/               Figure 6 (Membership Inference Attack)
    utility/               Figure 3 (Data Labelling) + supplementary GBMT
  chap/
    utility/               Figures 4-5 (TRTR/TSTR classifiers, GBMT clustering)
    variability/           Figure 7 (imputation vs. synthesis variance decomposition)
```

## Method → code map

| Paper element | Location |
|---|---|
| Algorithm 1 (cross-sectional eCDF-copula) | [`functions/synthesis_algorithms/algorithm1_cross_sectional.R`](functions/synthesis_algorithms/algorithm1_cross_sectional.R) |
| Algorithm 2 (longitudinal eCDF-copula) | [`functions/synthesis_algorithms/algorithm2_longitudinal.R`](functions/synthesis_algorithms/algorithm2_longitudinal.R) |
| REMBRANDT MICE imputation + long/wide reshape | [`functions/rembrandt/impute.R`](functions/rembrandt/impute.R) |
| REMBRANDT eCDF-copula synthesis | [`functions/rembrandt/generate_synthetic/eCDF_copula_generate.R`](functions/rembrandt/generate_synthetic/eCDF_copula_generate.R) |
| CHAP binning + MICE (M=5) | [`functions/chap/prepare_impute/`](functions/chap/prepare_impute/) |
| CHAP eCDF-copula + CTGAN ensemble (M=5 x L=10) | `02_mice_impute_and_ecdf_synthesize.R`, `CTGAN_generate.py` |
| Table 2 / Figures 1-2 (URA, MRA) | [`results/rembrandt/resemblance/`](results/rembrandt/resemblance/) |
| Figure 3 (Data Labelling Analysis) | [`results/rembrandt/resemblance/4_Data_Labelling_Resemblance_Dataset0.ipynb`](results/rembrandt/resemblance/4_Data_Labelling_Resemblance_Dataset0.ipynb) |
| Figure 4 (TRTR/TSTR, Dataset B) | [`results/chap/utility/`](results/chap/utility/) |
| Figure 5 (GBMT clustering ARI) | `results/chap/utility/gbmt_clustering/` |
| Figure 6 (Membership Inference Attack) | [`results/rembrandt/privacy/`](results/rembrandt/privacy/) |
| Figure 7 (variance decomposition) | [`results/chap/variability/`](results/chap/variability/) |

## Data availability

`data/raw/` contains the imputed REMBRANDT and CHAP patient data used to fit
the synthesizers. **This folder is excluded from version control (see
`.gitignore`) and is not part of the public repository** — it exists only in
the local copy of this project, since the underlying clinical data cannot be
shared publicly. Everything else (synthesis/imputation/evaluation code, the
synthetic datasets in `data/processed/`, and all evaluation results) is
public. Scripts that read from `data/raw/` will not run out of the box for
external users without access to the original REMBRANDT/CHAP data.

## Notes

- **Third-party code**: `functions/evaluation_functions/` vendors evaluation
  code from Hernandez et al. (2023) (MIT license). See
  `functions/evaluation_functions/NOTICE.md` for attribution.
- Notebooks use relative paths (`FUNCTIONS_HOME`, `REAL_DATA_HOME`,
  `SYN_DATA_HOME`) computed from each notebook's location — run them in place
  (don't move a notebook without updating these).
- GAN-based generation (CTGAN, WGAN-GP) requires Python (`ctgan`, `sdv`,
  `ydata-synthetic`); the eCDF-copula method and GBMT clustering require R
  (`copula`, `mice`, `gbmt`, `mclust`, `tidyverse`).
