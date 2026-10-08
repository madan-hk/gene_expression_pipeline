"""Checks run on the processed file BEFORE the file gets its real name."""
from __future__ import annotations

from pathlib import Path

import pandas as pd


class CheckError(RuntimeError):
    pass


def check_processed(path: Path, genes: list[str], sample_ids: list[str],
                    raw_sample_sums: pd.Series) -> dict:
    """Read the written file back and verify it. Raises CheckError on any problem."""
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")

    if df.index.name != "gene_id":
        raise CheckError("The first column is not called gene_id")
    if list(df.index) != genes:
        raise CheckError(f"Rows are not the {len(genes):,} protein-coding genes in annotation order")
    if list(df.columns) != sample_ids:
        raise CheckError("Sample columns or their order differ from the raw header")
    if df.isna().any().any():
        raise CheckError("There are empty values in the matrix")
    if (df.to_numpy() < 0).any():
        raise CheckError("There are negative values in the matrix")

    # Each sample's total here must equal the total of its kept genes in the raw file.
    # This catches a shifted, transposed or mixed-up matrix.
    got = df.sum(axis=0)
    want = raw_sample_sums.reindex(got.index)
    bad = got[got != want]
    if len(bad):
        raise CheckError(f"{len(bad)} samples have a different total than in the raw file (first: {bad.index[0]})")

    return {"n_samples": df.shape[1], "n_genes": df.shape[0], "dtype": str(df.dtypes.iloc[0])}
