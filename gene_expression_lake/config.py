"""Where things are, and the fixed choices (one annotation: G026; Human + bulk RNA-seq + blood)."""
from __future__ import annotations

import os
from pathlib import Path

ANNOTATION = "G026"          # GENCODE v26. Only this one is used.
RECOUNT3 = "https://recount-opendata.s3.amazonaws.com/recount3/release/human"

SOURCES = ("sra", "gtex", "tcga")

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# lake/ sits next to this file's project folder, whatever directory you run from.
# (RNA_LAKE_ROOT can point somewhere else, for example to try the pipeline in a scratch folder.)
LAKE_ROOT = Path(os.environ.get("RNA_LAKE_ROOT", PROJECT_ROOT / "lake"))

# The raw recount3 METADATA files (not the gene files): small files that describe each dataset's samples.
# They live inside the lake, in lake/metadata/<source>/<DATASET>/, and are used only to decide which datasets
# and samples are Human + bulk RNA-seq + blood. Metadata that sits in another folder is brought in with
# `add --from <folder>`, which copies it here.
METADATA_DIR = LAKE_ROOT / "metadata"

# The two lists of datasets this pipeline is allowed to build (made by the `select` command).
#   main list   : every sample of the dataset is whole blood and bulk
#   second list : every other dataset that has at least one bulk blood sample
SELECTION_DIR = Path(os.environ.get("RNA_SELECTION_DIR", PROJECT_ROOT / "selection"))
MAIN_LIST_FILE = SELECTION_DIR / "pure_whole_blood_datasets.tsv"
SECOND_LIST_FILE = SELECTION_DIR / "mixed_datasets.tsv"

# The 8 blood types and the folder of each. Whole blood has its own top-level folder; every
# other type goes under MIXED_BLOOD/.
BLOOD_TYPES = ["Whole blood", "PBMC", "Leukocytes", "Bone marrow", "Serum", "Plasma", "Buffy coat",
               "Blood (type not stated)"]
BLOOD_FOLDER = {
    "Whole blood": "WHOLE_BLOOD",
    "PBMC": "PBMC",
    "Leukocytes": "LEUKOCYTES",
    "Bone marrow": "BONE_MARROW",
    "Serum": "SERUM",
    "Plasma": "PLASMA",
    "Buffy coat": "BUFFY_COAT",
    "Blood (type not stated)": "BLOOD_TYPE_NOT_STATED",
}


def annotation_dir(root: Path = LAKE_ROOT) -> Path:
    return root / f"annotation_{ANNOTATION}"


def annotation_file(root: Path = LAKE_ROOT) -> Path:
    return annotation_dir(root) / f"human.gene_sums.{ANNOTATION}.gtf.gz"


def genes_file(root: Path = LAKE_ROOT) -> Path:
    """The protein-coding gene list, made once from the annotation file."""
    return annotation_dir(root) / "protein_coding_genes.tsv.gz"


def annotation_url() -> str:
    return f"{RECOUNT3}/annotations/gene_sums/human.gene_sums.{ANNOTATION}.gtf.gz"


def raw_dir(source: str, project: str, root: Path = LAKE_ROOT) -> Path:
    """raw/<PROJECT>/  (no source folder, same as processed/; the file name still starts with the source)"""
    return root / "raw" / project


def raw_file(source: str, project: str, root: Path = LAKE_ROOT) -> Path:
    return raw_dir(source, project, root) / f"{source}.gene_sums.{project}.{ANNOTATION}.gz"


def raw_url(source: str, project: str) -> str:
    # recount3 puts each project in a folder named by its last two characters.
    return (f"{RECOUNT3}/data_sources/{source}/gene_sums/{project[-2:]}/{project}/"
            f"{source}.gene_sums.{project}.{ANNOTATION}.gz")


def processed_dir(blood_type: str, source: str, project: str, root: Path = LAKE_ROOT) -> Path:
    """WHOLE_BLOOD/<PROJECT>/  or  MIXED_BLOOD/<BLOOD_TYPE>/<PROJECT>/

    There is no source folder (SRA / GTEX / TCGA): dataset ids are unique across the three sources
    (checked for all 818 datasets on the two lists), so the source is not needed to tell them apart.
    `source` stays in the signature because callers pass it.
    """
    folder = BLOOD_FOLDER[blood_type]
    if blood_type == "Whole blood":
        return root / "processed" / folder / project
    return root / "processed" / "MIXED_BLOOD" / folder / project


def processed_file(blood_type: str, source: str, project: str, root: Path = LAKE_ROOT) -> Path:
    return processed_dir(blood_type, source, project, root) / f"{project}.tsv.gz"
