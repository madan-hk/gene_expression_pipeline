# Gene Expression Lake

Builds the gene expression files of the RNA Annotation Workbench for **human, bulk RNA-seq, blood datasets** from
[recount3](https://rna.recount.bio/). For each qualifying dataset it keeps the **22,321 protein-coding genes**
(annotation G026, GENCODE v26) and only the bulk blood samples, and stores them in folders organised by blood type.

This repository holds the pipeline code and the dataset lists (`selection/`). The generated data (`lake/`) is
reproducible and excluded through `.gitignore`.

## Contents

1. [Setup](#setup)
2. [Quick start: generate the lake](#quick-start-generate-the-lake)
3. [Commands](#commands)
4. [Scope: the two dataset lists](#scope-the-two-dataset-lists)
5. [How a sample is judged](#how-a-sample-is-judged)
6. [How a dataset is built](#how-a-dataset-is-built)
7. [Output structure and file format](#output-structure-and-file-format)
8. [Storage and download size](#storage-and-download-size)
9. [Checking one dataset by hand](#checking-one-dataset-by-hand)
10. [Adding datasets](#adding-datasets)
11. [Safety and verification](#safety-and-verification)
12. [Limitations](#limitations)
13. [Troubleshooting](#troubleshooting)
14. [Project layout](#project-layout)

## Setup

Requirements: Python 3.10+, [uv](https://docs.astral.sh/uv/), internet access to recount3, and the recount3
**metadata** files (see below).

```bash
git clone https://github.com/madan-hk/gene_expression_pipeline.git && cd gene_expression_pipeline
uv sync
```

Run every command from the project folder.

**What it needs**

| Needed | Where it comes from |
|---|---|
| The code and `selection/` | the repository |
| The Python environment (`.venv`) | created by `uv sync` |
| The annotation and the 22,321-gene list | downloaded automatically on the first build (about 2 MB) |
| A dataset's raw gene file | downloaded automatically when that dataset is built |
| **The recount3 metadata** | **not downloaded by this project** (`lake/` is not in git, so a fresh clone has none). See below |

Without the metadata, `build` stops with a message naming the missing folder.

### Where the metadata is

The pipeline reads the metadata from `lake/metadata/`. Layout of a dataset's files:
`lake/metadata/<source>/<DATASET>/<source>.<kind>.<DATASET>.MD.gz` (for example
`lake/metadata/sra/SRP026331/sra.sra.SRP026331.MD.gz`).

**Metadata in another folder on your disk:**

- For the datasets **already on the lists** (for example after a fresh clone), copy the metadata folder in:
  `cp -R /path/to/metadata/. lake/metadata/`. `add` does not copy these, because it leaves listed datasets untouched.
- For **new** datasets, run `add --from /path/to/folder` (see [Adding datasets](#adding-datasets)). It searches that
  folder and its subfolders, puts the matching datasets on the main or second list, and **copies** their metadata into
  `lake/metadata/`.

## Quick start: generate the lake

From a fresh clone, run these in the project folder, in order. Every command is safe to run again: a finished dataset
is skipped, so an interrupted run continues where it stopped.

```bash
# 1. Install the environment
uv sync

# 2. Make the recount3 metadata available in lake/metadata/ (see "Where the metadata is" above):
#    - keep (or copy) your metadata folder there, with its sra/, gtex/, tcga/ folders:
#        cp -R /path/to/your/metadata/. lake/metadata/
#    - or, for datasets that are NOT on the lists yet, let `add` copy them in
#      (it creates lake/metadata/ if it does not exist):
#        uv run python -m gene_expression_lake.cli add --from /path/to/your/metadata

# 3. Check the two dataset lists
uv run python -m gene_expression_lake.cli list

# 4. Optional: see the download size before starting (slow)
uv run python -m gene_expression_lake.cli build-all --dry-run

# 5. Try one small dataset first
uv run python -m gene_expression_lake.cli build --source sra --project SRP158591

# 6. Build the main list (104 datasets), then the second list (718 datasets)
uv run python -m gene_expression_lake.cli build-all --list main
uv run python -m gene_expression_lake.cli build-all --list second
```

The first build also downloads the annotation and the gene list (about 2 MB). The results are in `lake/processed/`
(see [Output structure](#output-structure-and-file-format)) and the run logs in `lake/logs/`. None of this is tracked
by git.

## Commands

Prefix for every command: `uv run python -m gene_expression_lake.cli`

| What you want | Command |
|---|---|
| One dataset | `build --source sra --project SRP158591` |
| One dataset, only one blood type | `build --source sra --project SRP158591 --blood-type pbmc` |
| All datasets on the lists | `build-all` |
| Only the main list / only the second list | `build-all --list main` / `build-all --list second` |
| Only the first N datasets | `build-all --limit 10` |
| See the download size first | `build-all --dry-run` (slow: it asks the server for every file size) |
| Only some blood types | `build-all --blood-type pbmc` (repeatable) |
| Add new datasets from any folder | `add --from /path/to/folder` (`--dry-run` to preview) |
| Remake the lists from the metadata | `select` |
| Show the lists | `list` (`--top N` for more rows) |

`--source` is `sra`, `gtex` or `tcga`. A finished dataset is skipped on the next run (`--force` rebuilds it).
`build-all` builds 2 datasets at a time (`--workers N`), continues when one fails, and writes a log to `lake/logs/`.
Names for `--blood-type`: `whole_blood`, `pbmc`, `leukocytes`, `bone_marrow`, `serum`, `plasma`, `buffy_coat`,
`blood_type_not_stated`.

## Scope: the two dataset lists

Scope: **Human + bulk RNA-seq + any blood type**. A dataset qualifies when it is human **and has at least one sample
that is both bulk and blood**. After that it may have any number of other samples (single-cell, no prediction, not
blood), which are simply left out of the output.

Blood has 8 types: whole blood, PBMC, leukocytes, bone marrow, serum, plasma, buffy coat, and blood (type not stated).

| List | File | Rule |
|---|---|---|
| **Main** | `selection/pure_whole_blood_datasets.tsv` | **every sample** of the dataset is whole blood and bulk |
| **Second** | `selection/mixed_datasets.tsv` | every other qualifying dataset |

A dataset is in one list only. A dataset in neither list is refused, with the reason (for example *"no blood samples"*
or *"it has 40 blood samples but none are bulk RNA-seq"*). The files have a header line, so the number of datasets is
the number of lines minus one.

**Current state of the list files.** Made with `select` and the current rules (see below) on 8 Oct 2026: **822**
datasets out of the 8,867 recount3 projects of GTEx, TCGA and SRA: **104 main + 718 second** (3,777 samples in the
main list, 19,275 bulk blood samples in the second). Compared with the earlier lists (818 datasets, made with
`pred.type` only), the `curated.type` rule added 5 datasets (`SRP002105`, `SRP005279`, `SRP007338`, `SRP007569`,
`SRP050223`), removed 1 (`DRP001149`) and moved 5 from the second list to the main list (`SRP045500`, `SRP050000`,
`SRP056840`, `SRP059172`, `SRP062966`).

## How a sample is judged

The rules follow the workbench's sample filters (`backend/app/filters.py`).

| Question | Where it is read |
|---|---|
| Human | `organism` in `<source>.recount_project.<DATASET>.MD.gz` |
| Bulk RNA-seq | GTEx: always. TCGA: `gdc_experimental_strategy` is `RNA-Seq`. SRA, in `sra.recount_pred.<DATASET>.MD.gz`: **`curated.type` first** (set by hand, filled for only about 10% of samples); **when it is empty or `NA`, `pred.type`** (recount3's prediction). The value `rna-seq` = bulk, `scrna-seq` = single-cell, anything else = not bulk (reported as "no prediction") |
| Blood type | GTEx: `smtsd`. TCGA: `cgc_sample_sample_type` and `gdc_cases.project.primary_site`. SRA: tissue / source pairs in `sample_attributes`. Cell lines, blood vessels and culture-medium serum are not blood |

Notes on the SRA bulk rule, from the metadata of all 303,047 SRA samples:
- `curated.type` also holds values such as `target-locus`, `small-rna`, `miRNA`, `snRNA` and
  `hRNA-Seq(bacteria infected)`. All of these count as **not bulk**.
- `pred.type` and `pattern.predict.type` are different columns and often disagree. Only `curated.type` and `pred.type`
  are used here. `pred.type` says `scrna-seq` for many samples that are not single-cell (small RNA, target-locus).
- Bulk-versus-single-cell is a property of the sequencing. It is separate from the blood type (the source of the RNA).

## How a dataset is built

```text
 MAIN list     every sample is whole blood + bulk, so ALL samples are taken
                 -> WHOLE_BLOOD/<DATASET>/

 SECOND list   sample by sample:  keep bulk samples that have a blood type;
                                  leave out single-cell, no-prediction and no-blood samples
                 whole blood  -> WHOLE_BLOOD/<DATASET>/
                 other types  -> MIXED_BLOOD/<BLOOD_TYPE>/<DATASET>/

 neither       refused
```

Example, `SRP059205` (20 samples): 6 PBMC + 2 bone marrow + 3 blood-type-not-stated bulk samples are built (11); 2
single-cell, 2 with no prediction and 5 not-blood samples are left out.

## Output structure and file format

```text
lake/
├── metadata/<SOURCE>/<DATASET>/               recount3 metadata (input of `select`, `add` and `build`)
├── annotation_G026/                           the annotation and protein_coding_genes.tsv.gz (22,321 genes)
├── raw/<DATASET>/                             gene file exactly as recount3 provides it (name starts with the source)
├── processed/
│   ├── WHOLE_BLOOD/<DATASET>/<DATASET>.tsv.gz          (no source folder: dataset ids are unique across SRA, GTEx, TCGA)
│   └── MIXED_BLOOD/<BLOOD_TYPE>/<DATASET>/<DATASET>.tsv.gz
│                    (PBMC, LEUKOCYTES, BONE_MARROW, SERUM, PLASMA, BUFFY_COAT, BLOOD_TYPE_NOT_STATED)
└── logs/                                      one log per `build-all` run
```

A processed file is a compressed, tab-separated table: the first column is `gene_id`, every other column is one
sample, one row per gene (22,321 rows, same order in every file). Values are the raw recount3 coverage sums as whole
numbers (**not normalised**). Each file holds one blood type of one dataset. Open it with a double-click on macOS,
`gunzip -c <DATASET>.tsv.gz | head`, or `pandas.read_csv("<DATASET>.tsv.gz", sep="\t", index_col=0)`.

The source folder (SRA / GTEX / TCGA) is not part of the path. Dataset ids did not clash across the three sources in
any of the 822 datasets on the lists; this is not enforced by the code, so check it when adding datasets.

## Storage and download size

- Measured with `build-all --dry-run` on 8 Oct 2026, with the earlier 818-dataset lists: the raw gene files of the
  815 datasets not yet downloaded were **3.38 GB**. The 5 datasets added by the current rule are about 30 MB more, so
  expect **about 3.4 GB** for the 822 datasets. Processed files are smaller: only 22,321 of the genes and only the
  chosen samples.
- Little disk space? Build a part: `build-all --list main` (104 datasets), `build-all --limit 10`, or one dataset at a
  time.
- The files in `raw/` can be deleted to free space. A later build downloads them again.
- Very large datasets are read in slices, but have not been tested for memory use.

## Checking one dataset by hand

Take a small dataset and compare the output with the metadata.

1. Build it, for example `build --source sra --project SRP059205`.
2. Open the three metadata files in `lake/metadata/sra/<DATASET>/`:
   - `sra.recount_pred.<DATASET>.MD.gz`: `curated.type`, `pred.type` (bulk or single-cell).
   - `sra.sra.<DATASET>.MD.gz`: `external_id` (the sample id) and `sample_attributes` (tissue / source text).
   - `sra.recount_project.<DATASET>.MD.gz`: `organism`.
3. For each sample id decide: bulk (curated first, then `pred.type`) and blood type (from the tissue text).
4. The header of each output file should list exactly the samples you decided on, and no sample of another blood type.
5. The number of rows is 22,321, and the values are whole numbers.

## Adding datasets

Put the new datasets' metadata files **anywhere** on disk (any folder layout) and run:

```bash
uv run python -m gene_expression_lake.cli add --from /path/to/folder --dry-run   # preview
uv run python -m gene_expression_lake.cli add --from /path/to/folder             # add
uv run python -m gene_expression_lake.cli build --source sra --project <ID>      # build one of the new datasets
```

`build-all` builds **every** unfinished dataset on the lists, not only the new ones (about 815 datasets, 3.4 GB). To
build only what `add` just listed, use `build` for each new dataset.

`add` recognises a dataset by its sample file (`sra.sra.<ID>.MD.gz`, `gtex.gtex.<ID>.MD.gz`, `tcga.tcga.<ID>.MD.gz`)
and handles each one:

| Dataset | What `add` does |
|---|---|
| already on a list | left as it is (never added twice) |
| new and qualifying | metadata copied to `lake/metadata/<source>/<ID>/`, added to the right list |
| new, not qualifying | metadata copied, not added, reason printed |
| incomplete (needs `recount_project`; SRA also `recount_pred`) | not added, missing file named |

`select` rebuilds both lists from scratch from everything in the metadata folder (use it after changing the rules or
the metadata). It overwrites the list files, so back them up first (`cp -R selection selection_backup`; ignored by git).
The lists are generated: do not edit them by hand.
If a dataset moves lists or drops out, its earlier files stay in `lake/processed/` and must be deleted by hand.
A new source (not GTEx, TCGA or SRA) needs its rules added in `gene_expression_lake/selection.py`.

## Safety and verification

- Downloads and output files are written under a temporary `.part` name and renamed only when complete; output files
  only after the checks pass, so a file with its real name is always complete.
- Every output file is checked: exactly the 22,321 protein-coding genes in the fixed order, exactly the chosen
  samples, no empty or negative values, and each sample's total equals its total in the raw file.
- The sample choice (per dataset and blood type) was compared with the workbench's ingestion adapters for all 818
  datasets of the earlier lists with no differences. That comparison was made **before** the `curated.type` rule was
  added, so it has not been repeated for the current rule.
- The current list files are the output of `select` over the metadata of all 8,867 projects (822 datasets).

## Limitations

- **iPSC-derived samples can count as blood.** The blood type is read from keywords such as "PBMC" or "cord blood" in the
  `source_name` text. Lab-grown induced pluripotent stem cell (iPSC) lines made *from* PBMC or cord blood carry those
  words, so they are classed as blood even though the sequenced cells are not. In `SRP059205`, 9 of the 11 samples
  built are of this kind. This is not corrected yet.
- Blood type is derived from the words in a sample's tissue text, so unusual wording can be missed.
- SRA bulk-versus-single-cell comes from recount3's hand-set `curated.type` and prediction `pred.type`, and can be wrong.
- The values are raw coverage sums. Normalisation is not part of this pipeline.
- Only the protein-coding genes of annotation G026 are kept. It has not been compared with another annotation (G029).
- Dataset ids are assumed unique across SRA, GTEx and TCGA (see above).

## Troubleshooting

Each problem below shows what you see, why it happens, and what to do. Run all commands from the project folder.

**1. `FAILED SRA/<ID>: its metadata folder was not found at ...`**
- Why: the metadata files of that dataset are not in `lake/metadata/<source>/<ID>/`. The pipeline needs them to decide
  which samples are bulk and blood.
- Do: put that dataset's metadata files in the folder shown in the message, or run
  `add --from /folder/that/has/them`, which copies them there.

**2. `FAILED ... is not built: <reason>. (It is in neither list.)`**
- Why: the dataset is on neither the main list nor the second list, so it is not allowed to be built. The reason is
  printed, for example "not human", "no blood samples", or "it has 40 blood samples but none are bulk RNA-seq".
- Do: if the reason is right, nothing. The dataset does not qualify. If it is a new dataset, add it with
  `add --from <folder>`. If you changed the rules or the metadata since the lists were made, remake the lists with
  `select` (back up the `selection/` folder first, because `select` overwrites both files).

**3. `FAILED ... is in the main list but its metadata no longer says every sample is whole blood + bulk. Run select again.`**
- Why: the list file says every sample of this dataset is whole blood and bulk, but the metadata, read with the current
  rules, says otherwise. This happens after the rules or the metadata changed and the lists were not remade.
- Do: back up `selection/`, then run `select`. The dataset moves to the second list or drops out.

**4. `Metadata folder not found: ...` when running `select`**
- Why: `lake/metadata/` does not exist, so there is nothing to scan.
- Do: copy the recount3 metadata into `lake/metadata/`, or run `add --from <folder>` for the datasets you have.

**5. `Not found on recount3: <url>`**
- Why: recount3 has no gene file at that address for this dataset (G026 annotation). The dataset id may be wrong or
  missing from recount3.
- Do: check the dataset id and source (`--source sra`, `gtex` or `tcga`).

**6. `Download failed after 3 attempts`**
- Why: network problem, or recount3 is not reachable right now.
- Do: run the same command again. Finished files are kept and skipped, and a half-downloaded file is thrown away.

**7. A file ending in `.part` inside `lake/`**
- Why: a run was stopped (or crashed) while writing that file. A real file name only appears when the file is complete.
- Do: nothing. The next run replaces it. You can also delete it.

**8. `warning: VIRTUAL_ENV=... does not match the project environment`**
- Why: another Python environment is active in your terminal. `uv` ignores it and uses the project's own.
- Do: nothing, it is harmless. To silence it, run `deactivate` or start commands with `env -u VIRTUAL_ENV`.

**9. `No module named gene_expression_lake`, or `uv: command not found`**
- Why: you are not in the project folder, the packages are not installed yet, or `uv` is not installed.
- Do: `cd` to the project folder and run `uv sync` once. Install `uv` first if the command is missing.

**10. `build-all --dry-run` takes many minutes**
- Why: it asks recount3 for the size of each file, one by one (about 815 files).
- Do: wait, or use `build-all --list main --dry-run` or `--limit 10` to check fewer files.

## Project layout

```text
gene_expression_lake/
├── config.py       paths, settings (RNA_* variables), annotation, list files, blood-type folders
├── selection.py    sorts samples by blood and sequencing type; makes, updates and reads the lists
├── annotation.py   the protein-coding gene list
├── download.py     safe downloads
├── builder.py      builds the genes x samples table from a raw gene file
├── checks.py       verifies an output file before it is accepted
└── cli.py          commands: add, select, list, build, build-all
selection/          the two dataset lists
```
