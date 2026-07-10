# =============================================================================
# CHAP: CTGAN Synthesis Pipeline
# Run AFTER: 2_mice_ecdf.R  (which saves output/chap/train/train_imp_0X.csv)
#
# Output folder structure:
#   output/chap/
#     synthetic_ctgan/  imp_01/syn_01.csv … imp_05/syn_10.csv
#     metadata/         run_log_ctgan.csv
#
# Seed strategy (fully deterministic from 3 constants below):
#   ctgan fit   : SEED_CTGAN_BASE  + imp_id          (e.g. 3001 … 3005)
#   ctgan sample: SEED_SAMPLE_BASE + (imp_id-1)*N_SYN + syn_id  (e.g. 4001 … 4050)
# =============================================================================

import warnings
warnings.filterwarnings("ignore")

import gc
import os
import numpy as np
import pandas as pd
import torch
from ctgan import CTGAN
from datetime import datetime

# ── 1. Constants ──────────────────────────────────────────────────────────────
SEED_CTGAN_BASE  = 3000   # fit seed for imp i  →  3000 + i
SEED_SAMPLE_BASE = 4000   # sample seed for (imp i, syn j) → 4000 + (i-1)*N_SYN + j

N_IMP = 5
N_SYN = 10

BASE_OUT  = "./output/chap"
TRAIN_DIR = os.path.join(BASE_OUT, "train")
SYN_DIR   = os.path.join(BASE_OUT, "synthetic_ctgan")
META_DIR  = os.path.join(BASE_OUT, "metadata")

DISCRETE_COLUMNS = ["Group"]

# ── 2. Setup Output Directories ───────────────────────────────────────────────
for i in range(1, N_IMP + 1):
    os.makedirs(os.path.join(SYN_DIR, f"imp_{i:02d}"), exist_ok=True)
os.makedirs(META_DIR, exist_ok=True)

# ── 3. Run Pipeline ───────────────────────────────────────────────────────────
log_path = os.path.join(META_DIR, "run_log_ctgan.csv")
log_header_written = os.path.exists(log_path)

for i in range(1, N_IMP + 1):
    imp_tag    = f"imp_{i:02d}"
    train_file = os.path.join(TRAIN_DIR, f"train_{imp_tag}.csv")
    fit_seed   = SEED_CTGAN_BASE + i

    print(f"\n══ Imputation {i} (fit seed {fit_seed}) ════════════════════════")

    train_df = pd.read_csv(train_file)
    n_train  = len(train_df)
    print(f"  Loaded {train_file}: {n_train} rows")

    # Fit CTGAN (set global seeds for reproducibility)
    np.random.seed(fit_seed)
    torch.manual_seed(fit_seed)
    ctgan = CTGAN()
    ctgan.fit(train_df, DISCRETE_COLUMNS)
    del train_df
    gc.collect()
    print(f"  CTGAN fit done.")

    # Sample → save → free, one dataset at a time
    for j in range(1, N_SYN + 1):
        sample_seed = SEED_SAMPLE_BASE + (i - 1) * N_SYN + j
        np.random.seed(sample_seed)
        torch.manual_seed(sample_seed)

        syn_df   = ctgan.sample(n_train)
        syn_file = os.path.join(SYN_DIR, imp_tag, f"syn_{j:02d}.csv")
        syn_df.to_csv(syn_file, index=False)

        # Append one row to log immediately, no accumulation in memory
        log_row = pd.DataFrame([{
            "stage":       "ctgan_syn",
            "imp_id":      i,
            "syn_id":      j,
            "fit_seed":    fit_seed,
            "sample_seed": sample_seed,
            "n_rows":      len(syn_df),
            "n_group0":    (syn_df["Group"] == "Group0").sum(),
            "n_group1":    (syn_df["Group"] == "Group1").sum(),
            "file_path":   os.path.abspath(syn_file),
            "timestamp":   datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }])
        log_row.to_csv(log_path, mode="a", index=False, header=not log_header_written)
        log_header_written = True

        del syn_df, log_row
        gc.collect()

        print(f"  syn {j:02d} (sample seed {sample_seed}) → {os.path.basename(syn_file)}")

    # Free the model before loading the next imputation
    del ctgan
    gc.collect()

print(f"\n══ Done ══════════════════════════════════════════════════════")
print(f"Imputed datasets  : {N_IMP}")
print(f"Synthetic datasets: {N_IMP * N_SYN}  ({N_IMP} imp × {N_SYN} syn)")
print(f"run_log saved → {log_path}")
