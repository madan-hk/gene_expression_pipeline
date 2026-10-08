"""Turn one raw recount3 gene file into the processed file.

raw:        63,856 genes as rows, samples as columns, with ## comment lines
processed:  the 22,321 protein-coding genes as rows (annotation order) and ONLY the chosen samples
            as columns, as a tab-separated text file (gzip-compressed): gene_id, then one column per sample
"""
from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

CHUNK_ROWS = 10_000   # read the raw file in slices, so the whole 63,856-gene table is never in memory


class BuildError(RuntimeError):
    pass


@dataclass
class Built:
    matrix: pd.DataFrame          # index = gene_id, one column per sample
    sample_ids: list[str]         # the samples in the raw header, in order
    raw_sample_sums: pd.Series    # per-sample total of the kept genes, counted while reading the raw file

    def subset(self, ids: list[str]) -> "Built":
        """The same table with only these samples (columns keep the raw file's order)."""
        keep = set(ids)
        cols = [c for c in self.sample_ids if c in keep]
        return Built(matrix=self.matrix[cols], sample_ids=cols, raw_sample_sums=self.raw_sample_sums.reindex(cols))


def _header_samples(raw_path: Path) -> list[str]:
    """Sample ids from the raw header line (the first line that is not a ## comment)."""
    with gzip.open(raw_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if not cols or cols[0] != "gene_id":
                raise BuildError(f"Unexpected header in {raw_path.name}: starts with {cols[0]!r}, expected 'gene_id'")
            return cols[1:]
    raise BuildError(f"{raw_path.name} has no header line")


def build_matrix(raw_path: Path, genes: list[str], chosen: list[str] | None = None) -> Built:
    """`chosen`: the sample ids to keep, or None for every sample of the raw file.
    The columns come out in the raw file's own order."""
    header = _header_samples(raw_path)
    if len(set(header)) != len(header):
        raise BuildError("The raw file has duplicate sample ids in its header")
    if chosen is None:
        samples = list(header)
    else:
        wanted_samples = set(chosen)
        absent = wanted_samples - set(header)
        if absent:
            raise BuildError(f"{len(absent)} chosen samples are not in the raw gene file (first: {sorted(absent)[0]})")
        samples = [h for h in header if h in wanted_samples]
    if not samples:
        raise BuildError("No chosen samples to build")

    wanted = set(genes)
    kept: list[pd.DataFrame] = []
    sums = pd.Series(0, index=samples, dtype="int64")
    for chunk in pd.read_csv(raw_path, sep="\t", comment="#", index_col=0, usecols=["gene_id"] + samples, chunksize=CHUNK_ROWS):
        part = chunk[chunk.index.isin(wanted)]
        if len(part):
            kept.append(part)
            sums = sums.add(part.sum(axis=0), fill_value=0)
    if not kept:
        raise BuildError("None of the protein-coding genes were found in the raw file")

    table = pd.concat(kept)                     # kept genes x samples
    if not table.index.is_unique:
        raise BuildError("The raw file lists some gene more than once")
    missing = [g for g in genes if g not in table.index]
    if missing:
        raise BuildError(f"{len(missing)} protein-coding genes are missing from the raw file (first: {missing[0]})")

    table = table.loc[genes]                    # annotation order, the same for every dataset
    table.index.name = "gene_id"
    return Built(matrix=table, sample_ids=samples, raw_sample_sums=sums.astype("int64"))


def write_tsv_gz(matrix: pd.DataFrame, path: Path) -> int:
    """Write the matrix as gzip-compressed TSV at `path` (the caller gives a temporary name).

    The name stored inside the gzip is `expression.tsv`, so unzipping it on a Mac or with
    gunzip gives a file called expression.tsv, not the temporary name. Returns the size in bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="expression.tsv", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as text:
                matrix.to_csv(text, sep="\t")
    return path.stat().st_size
