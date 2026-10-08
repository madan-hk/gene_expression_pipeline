"""The G026 annotation: one file for the whole lake, and the protein-coding gene list made from it.

The list is made ONCE and saved as lake/annotation_G026/protein_coding_genes.tsv.gz
(columns: gene_id, gene_name). Every later run reads that small file; the big annotation
file is only needed (and downloaded) the first time, to make it.
"""
from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

from . import config
from .download import download

_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_GENE_TYPE = re.compile(r'gene_type "([^"]+)"')
_GENE_NAME = re.compile(r'gene_name "([^"]+)"')
_HEADER = "gene_id\tgene_name"


def ensure_annotation(root: Path = config.LAKE_ROOT) -> Path:
    """Download the annotation file once. Later calls find it and do nothing."""
    path = config.annotation_file(root)
    download(config.annotation_url(), path)
    return path


def protein_coding_genes(path: Path) -> list[tuple[str, str]]:
    """(gene_id, gene_name) of every gene whose gene_type is protein_coding, in the order of
    the annotation file.

    The order is fixed by the file, so every dataset's processed matrix has the same
    genes (rows) in the same order."""
    genes: list[tuple[str, str]] = []
    seen: set[str] = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            t = _GENE_TYPE.search(cols[8])
            if not t or t.group(1) != "protein_coding":
                continue
            g = _GENE_ID.search(cols[8])
            if g and g.group(1) not in seen:
                seen.add(g.group(1))
                n = _GENE_NAME.search(cols[8])
                genes.append((g.group(1), n.group(1) if n else ""))
    if not genes:
        raise ValueError(f"No protein_coding genes found in {path}")
    return genes


def _write_genes(genes: list[tuple[str, str]], target: Path) -> None:
    """Save the list as gzip-compressed TSV, under a temporary name first. The name stored inside
    the gzip is protein_coding_genes.tsv, so double-clicking it on a Mac unpacks to that name."""
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    try:
        with open(part, "wb") as raw:
            with gzip.GzipFile(filename="protein_coding_genes.tsv", mode="wb", fileobj=raw,
                               compresslevel=6, mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as text:
                    text.write(_HEADER + "\n")
                    for gene_id, gene_name in genes:
                        text.write(f"{gene_id}\t{gene_name}\n")
        part.replace(target)
    finally:
        part.unlink(missing_ok=True)


def _read_genes(path: Path) -> list[str] | None:
    """The gene ids from the saved list, or None if the file is missing or damaged."""
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt") as f:
            if f.readline().rstrip("\n") != _HEADER:
                return None
            ids = [line.split("\t", 1)[0] for line in f if line.strip()]
    except (OSError, EOFError):
        return None
    if not ids or len(set(ids)) != len(ids):
        return None
    return ids


def load_genes(root: Path = config.LAKE_ROOT) -> list[str]:
    """The protein-coding gene ids, in the fixed order.

    Reads the saved list if it is there and valid. Otherwise downloads the annotation (if it is not
    already there), makes the list, saves it, and returns it."""
    saved = config.genes_file(root)
    ids = _read_genes(saved)
    if ids is not None:
        return ids
    genes = protein_coding_genes(ensure_annotation(root))
    _write_genes(genes, saved)
    return [gene_id for gene_id, _ in genes]
