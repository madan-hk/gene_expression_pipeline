"""Which datasets and samples are in scope: Human + bulk RNA-seq + blood (any of 8 blood types).

The rules are the same as the workbench's sample filters (backend/app/filters.py):
  * Human          organism is Homo sapiens (the recount_project file)
  * Bulk RNA-seq   GTEx: always.  TCGA: gdc_experimental_strategy is "RNA-Seq".
                   SRA: curated.type when it has a value (set by hand), otherwise recount3's pred.type;
                   "rna-seq" = bulk, "scrna-seq" = single-cell, anything else = no prediction
  * Blood type     the sample's own tissue words, sorted into 8 types in a fixed order (PBMC, buffy coat,
                   whole blood, plasma, serum, leukocytes, bone marrow, blood type not stated). Lab cell lines,
                   blood vessels and culture-medium serum are not blood.

Reads only the raw METADATA files (config.METADATA_DIR), never the gene files.

Two lists of datasets (made once by `select`, then only read):
  * main list    every sample of the dataset is whole blood and bulk       -> all samples are taken directly
  * second list  every other dataset with at least one bulk blood sample   -> sample by sample
A dataset in neither list is refused.
"""
from __future__ import annotations

import csv
import os
import re
import shutil
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import config

# Blood types in the same priority order as the workbench; the first match wins.
_BLOOD = [
    ("PBMC", r"\b(pbmcs?|peripheral blood mononuclear)"),
    ("Buffy coat", r"\bbuffy"),
    ("Whole blood", r"\b(whole blood|peripheral blood)\b"),
    ("Plasma", r"\bplasma\b"),
    ("Serum", r"\bserum\b"),
    ("Leukocytes", r"\b(leukocytes?|lymphocytes?|monocytes?|neutrophils?|granulocytes?|t cells?|b cells?|nk cells?|white blood cells?)\b"),
    ("Bone marrow", r"\bbone marrow\b"),
    ("Blood (type not stated)", r"\bblood\b"),
]
_BLOOD = [(name, re.compile(rx)) for name, rx in _BLOOD]
_NOT_BLOOD = re.compile(r"transformed|cell line|blood vessel|blood[- ]brain")
_CULTURE = re.compile(r"medium|media|cultur|grown|bovine|fbs|fetal calf")
_ATTR_TISSUE = re.compile(r"(?:^|\|)(?:tissue|source_name|specimen|sample[ _]type|body[ _]site);;([^|]*)", re.I)
_HINT_KEYS = ["tissue", "source_name", "cell_type", "cell_line", "biomaterial"]

BULK, SINGLE, NO_PRED = "bulk", "single-cell", "no prediction"


def blood_type(text: str) -> str | None:
    t = text.lower()
    if _NOT_BLOOD.search(t):
        return None
    for name, rx in _BLOOD:
        if rx.search(t):
            if name in ("Plasma", "Serum") and _CULTURE.search(t):
                continue
            return name
    return None


def _tissue_hint(attributes: str) -> str:
    """First value in SRA's packed attributes whose key looks like a tissue / source field."""
    pairs: dict[str, str] = {}
    for chunk in str(attributes).split("|"):
        if ";;" in chunk:
            k, v = chunk.split(";;", 1)
            pairs[k.strip().lower()] = v.strip()
    for kw in _HINT_KEYS:
        for k, v in pairs.items():
            if kw in k:
                return v
    return ""


def _sra_sequencing(curated: str | None, predicted: str | None) -> str:
    """Bulk, single-cell or neither for one SRA sample, from two columns of the recount_pred file:
      1. `curated.type`  filled by hand for only about 10% of samples; used first when it has a value.
                         "rna-seq" = bulk, "scrna-seq" = single-cell, any other value (target-locus, small-rna,
                         miRNA, snRNA, ...) = not bulk, reported as NO_PRED.
      2. `pred.type`     recount3's prediction, used only when `curated.type` is empty or "NA".
                         "rna-seq" = bulk, "scrna-seq" = single-cell, anything else = NO_PRED.
    """
    value = curated if curated not in (None, "", "NA") else predicted
    return BULK if value == "rna-seq" else SINGLE if value == "scrna-seq" else NO_PRED


def _organism(folder: Path, source: str, project: str) -> str:
    f = folder / f"{source}.recount_project.{project}.MD.gz"
    if not f.exists():
        return ""
    d = pd.read_csv(f, sep="\t", dtype=str, nrows=1, keep_default_na=False)
    return d["organism"].iloc[0] if "organism" in d else ""


@dataclass
class Summary:
    """What one dataset's metadata says. `samples` is (sample id, blood type or None, sequencing)."""
    source: str
    project: str
    human: bool
    samples: list[tuple[str, str | None, str]] = field(default_factory=list)

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    def count(self, blood: str | None = "any", seq: str | None = None) -> int:
        """Samples with a blood type ("any" = any type, None = not blood, or a type name) and a sequencing."""
        return sum(1 for _, b, s in self.samples
                   if (blood == "any" and b is not None or b == blood and blood != "any")
                   and (seq is None or s == seq))

    @property
    def n_bulk_blood(self) -> int:
        return self.count("any", BULK)

    @property
    def n_not_blood(self) -> int:
        return sum(1 for _, b, _ in self.samples if b is None)

    def bulk_by_type(self) -> dict[str, list[str]]:
        """{blood type: [sample ids]} for the bulk blood samples, in the order of the sample file."""
        out: dict[str, list[str]] = {}
        for sid, b, s in self.samples:
            if b is not None and s == BULK:
                out.setdefault(b, []).append(sid)
        return out

    @property
    def is_main(self) -> bool:
        """Every sample of the dataset is whole blood and bulk."""
        return self.human and self.n_samples > 0 and self.count("Whole blood", BULK) == self.n_samples

    @property
    def in_scope(self) -> bool:
        return self.human and self.n_bulk_blood > 0

    @property
    def group(self) -> str:
        """Kind of dataset, judged on its blood samples."""
        single, nopred = self.count("any", SINGLE) > 0, self.count("any", NO_PRED) > 0
        if not single and not nopred:
            return "pure bulk"
        return "bulk + " + " + ".join(x for x, on in (("single-cell", single), ("no prediction", nopred)) if on)


def summarize(source: str, project: str, metadata_dir: Path | None = None) -> Summary:
    """Read one dataset's metadata and sort every sample (blood type, sequencing)."""
    folder = (metadata_dir or config.METADATA_DIR) / source / project
    human = _organism(folder, source, project).strip().lower() == "homo sapiens"
    sm = Summary(source, project, human)
    if not human:
        return sm

    if source == "sra":
        sf = folder / f"sra.sra.{project}.MD.gz"
        if not sf.exists():
            return sm
        df = pd.read_csv(sf, sep="\t", dtype=str, keep_default_na=False,
                         usecols=lambda c: c.lower() in ("external_id", "sample_attributes"))
        df.columns = [c.lower() for c in df.columns]
        pf = folder / f"sra.recount_pred.{project}.MD.gz"
        pred: dict[str, str] = {}
        curated: dict[str, str] = {}
        if pf.exists():
            pr = pd.read_csv(pf, sep="\t", dtype=str, keep_default_na=False,
                             usecols=lambda c: c in ("external_id", "pred.type", "curated.type"))
            pred = dict(zip(pr["external_id"], pr["pred.type"]))
            if "curated.type" in pr:
                curated = dict(zip(pr["external_id"], pr["curated.type"]))
        attrs = df["sample_attributes"] if "sample_attributes" in df else [""] * len(df)
        for ext, a in zip(df["external_id"], attrs):
            m = _ATTR_TISSUE.search(a)
            bt = blood_type(_tissue_hint(a) + " " + (m.group(1) if m else ""))
            sm.samples.append((ext, bt, _sra_sequencing(curated.get(ext), pred.get(ext))))
        return sm

    sf = folder / f"{source}.{source}.{project}.MD.gz"
    if not sf.exists():
        return sm
    df = pd.read_csv(sf, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
    df.columns = [c.lower() for c in df.columns]
    for _, r in df.iterrows():
        if source == "gtex":
            text, seq = r.get("smtsd", ""), BULK
        else:
            text = r.get("cgc_sample_sample_type", "") + " " + r.get("gdc_cases.project.primary_site", "")
            seq = BULK if r.get("gdc_experimental_strategy", "") == "RNA-Seq" else NO_PRED
        sm.samples.append((r["external_id"], blood_type(text), seq))
    return sm


def refusal_reason(sm: Summary) -> str:
    """Why a dataset is in neither list, in plain words."""
    if sm.n_samples == 0 and not sm.human:
        return "not found in the metadata, or not human"
    if not sm.human:
        return "not human"
    if sm.n_samples == 0:
        return "no samples found in the metadata"
    blood = sm.count("any")
    if blood == 0:
        return f"no blood samples (0 of {sm.n_samples} samples have a blood type)"
    return (f"it has {blood} blood samples but none are bulk RNA-seq "
            f"({sm.count('any', SINGLE)} single-cell, {sm.count('any', NO_PRED)} with no prediction)")


def _scan_one(job: tuple[str, str]) -> tuple[str, str, Summary | None]:
    source, project = job
    try:
        return source, project, summarize(source, project)
    except Exception:
        return source, project, None


MAIN_COLUMNS = ["source", "project", "n_samples"]
SECOND_COLUMNS = (["source", "project", "group", "n_bulk_blood", "n_single_cell_blood",
                   "n_no_prediction_blood", "n_not_blood"]
                  + ["bulk_" + b.lower().replace(" ", "_").replace("(", "").replace(")", "") for b in config.BLOOD_TYPES])


def _write(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def row_for(sm: Summary) -> tuple[str, dict] | None:
    """("main" or "second", the row for that list) for a dataset in scope, or None if it is not in scope."""
    if sm is None or not sm.in_scope:
        return None
    if sm.is_main:
        return "main", {"source": sm.source.upper(), "project": sm.project, "n_samples": sm.n_samples}
    row = {"source": sm.source.upper(), "project": sm.project, "group": sm.group, "n_bulk_blood": sm.n_bulk_blood,
           "n_single_cell_blood": sm.count("any", SINGLE), "n_no_prediction_blood": sm.count("any", NO_PRED),
           "n_not_blood": sm.n_not_blood}
    for b, col in zip(config.BLOOD_TYPES, SECOND_COLUMNS[7:]):
        row[col] = sm.count(b, BULK)
    return "second", row


def _sorted(main_rows: list[dict], second_rows: list[dict]) -> tuple[list[dict], list[dict]]:
    main_rows = sorted(main_rows, key=lambda r: (-int(r["n_samples"]), r["source"], r["project"]))
    second_rows = sorted(second_rows, key=lambda r: (-int(r["n_bulk_blood"]), r["source"], r["project"]))
    return main_rows, second_rows


def build_lists(workers: int | None = None) -> tuple[list[dict], list[dict]]:
    """Scan every project in the metadata folder and write the two lists from scratch. Returns (main, second)."""
    md = config.METADATA_DIR
    if not md.exists():
        raise SystemExit(f"Metadata folder not found: {md}\nPut the recount3 metadata in lake/metadata/ (datasets in another folder can be brought in with `add --from <folder>`).")
    jobs = [(s, p.name) for s in config.SOURCES for p in sorted((md / s).glob("*")) if p.is_dir()]
    with ProcessPoolExecutor(max_workers=workers or os.cpu_count()) as ex:
        results = list(ex.map(_scan_one, jobs, chunksize=40))
    failed = [(s, p) for s, p, sm in results if sm is None]
    main_rows, second_rows = [], []
    for s, p, sm in results:
        placed = row_for(sm)
        if placed:
            (main_rows if placed[0] == "main" else second_rows).append(placed[1])
    main_rows, second_rows = _sorted(main_rows, second_rows)
    _write(config.MAIN_LIST_FILE, MAIN_COLUMNS, main_rows)
    _write(config.SECOND_LIST_FILE, SECOND_COLUMNS, second_rows)
    if failed:
        print(f"Warning: {len(failed)} projects could not be read, for example {failed[:3]}")
    return main_rows, second_rows


# ---------------------------------------------------------------------------------------------
# Adding datasets from any folder, without touching the datasets already on the lists
# ---------------------------------------------------------------------------------------------
_SAMPLE_FILE = re.compile(r"^(sra|gtex|tcga)\.(sra|gtex|tcga)\.([A-Za-z0-9_]+)\.MD\.gz$")
_META_KINDS = ("recount_project", "recount_pred", "recount_qc", "recount_seq_qc")


def find_datasets(folder: Path) -> dict[tuple[str, str], list[Path]]:
    """Find the datasets inside ANY folder, whatever its layout: every `<source>.<source>.<DATASET>.MD.gz` file
    is a dataset. Returns {(source, project): [its metadata files found next to it]}."""
    found: dict[tuple[str, str], list[Path]] = {}
    for f in sorted(folder.rglob("*.MD.gz")):
        m = _SAMPLE_FILE.match(f.name)
        if not m or m.group(1) != m.group(2):
            continue
        source, project = m.group(1), m.group(3)
        files = [f] + [f.with_name(f"{source}.{k}.{project}.MD.gz") for k in _META_KINDS
                       if f.with_name(f"{source}.{k}.{project}.MD.gz").exists()]
        found.setdefault((source, project), files)       # the first copy found wins
    return found


def copy_into_metadata(source: str, project: str, files: list[Path]) -> list[str]:
    """Copy a dataset's metadata files into lake/metadata/<source>/<DATASET>/ (never overwriting). Returns the
    names of the files copied."""
    dest = config.METADATA_DIR / source / project
    dest.mkdir(parents=True, exist_ok=True)
    copied = []
    for f in files:
        target = dest / f.name
        if not target.exists():
            shutil.copy2(f, target)
            copied.append(f.name)
    return copied


def read_rows(path: Path) -> list[dict]:
    return _read_list(path)


def save_lists(main_rows: list[dict], second_rows: list[dict]) -> None:
    main_rows, second_rows = _sorted(main_rows, second_rows)
    _write(config.MAIN_LIST_FILE, MAIN_COLUMNS, main_rows)
    _write(config.SECOND_LIST_FILE, SECOND_COLUMNS, second_rows)


def _read_list(path: Path) -> list[dict]:
    if not path.exists():
        raise SystemExit("The lists of datasets are missing. Run:  uv run python -m gene_expression_lake.cli select")
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_lists() -> dict[tuple[str, str], str]:
    """{(source, project): "main" or "second"} from the saved lists. A dataset in neither is not in the dict."""
    out: dict[tuple[str, str], str] = {}
    for r in _read_list(config.MAIN_LIST_FILE):
        out[(r["source"].lower(), r["project"])] = "main"
    for r in _read_list(config.SECOND_LIST_FILE):
        out[(r["source"].lower(), r["project"])] = "second"
    return out
