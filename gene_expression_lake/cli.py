"""Build the gene expression files for Human + bulk RNA-seq + blood datasets.

    uv run python -m gene_expression_lake.cli add --from /any/folder       # add new datasets to the lists (from any folder)
    uv run python -m gene_expression_lake.cli select                       # make the two lists again from lake/metadata/
    uv run python -m gene_expression_lake.cli list                         # show them
    uv run python -m gene_expression_lake.cli build --source sra --project SRP158591     # ONE dataset
    uv run python -m gene_expression_lake.cli build-all --dry-run          # how much would be downloaded
    uv run python -m gene_expression_lake.cli build-all                    # every dataset in the lists

How a dataset is built depends on the list it is in:
  main list    every sample is whole blood + bulk, so ALL samples are taken -> WHOLE_BLOOD/<SOURCE>/<DATASET>/
  second list  sample by sample: only bulk samples that have a blood type are kept (single-cell, no
               prediction and no-blood samples are left out), split by blood type:
                 whole blood -> WHOLE_BLOOD/<SOURCE>/<DATASET>/
                 other types -> MIXED_BLOOD/<BLOOD_TYPE>/<SOURCE>/<DATASET>/
  neither      refused, with the reason

Every file is written under a temporary name, checked, and only then renamed to <DATASET>.tsv.gz, so a file
with that name is always complete. A finished file is skipped unless --force.
"""
from __future__ import annotations

import argparse
import re
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path

from . import config
from .annotation import load_genes
from .builder import BuildError, build_matrix, write_tsv_gz
from .checks import CheckError, check_processed
from .download import DownloadError, download, remote_size
from .selection import (build_lists, copy_into_metadata, find_datasets, load_lists, read_rows, refusal_reason,
                        row_for, save_lists, summarize)

_TYPE_NAMES = {folder.lower(): bt for bt, folder in config.BLOOD_FOLDER.items()}   # "pbmc" -> "PBMC", ...


def _say(step: str, text: str) -> None:
    print(f"[{step:<10}] {text}", flush=True)


def build_one(source: str, project: str, force: bool = False, quiet: bool = False,
              only_types: set[str] | None = None) -> tuple[str, str]:
    """Build one dataset. Returns (status, message); status is done, skipped or failed."""
    say = (lambda *a: None) if quiet else _say
    source = source.lower()
    if source not in config.SOURCES:
        return "failed", f"Unknown source {source!r}. Use one of: {', '.join(config.SOURCES)}"
    if not re.fullmatch(r"[A-Za-z0-9_]+", project):
        return "failed", f"Project name {project!r} may only contain letters, digits and underscores"
    name = f"{source.upper()}/{project}"

    if not (config.METADATA_DIR / source / project).is_dir():
        return "failed", (f"{name}: its metadata folder was not found at {config.METADATA_DIR / source / project}. "
                          f"Put the recount3 metadata in {config.METADATA_DIR}, or bring it in with `add --from <folder>`.")

    kind = load_lists().get((source, project))
    sm = summarize(source, project)
    if kind is None:
        return "failed", f"{name} is not built: {refusal_reason(sm)}. (It is in neither list.)"

    # What to build: {blood type: sample ids}; None means every sample of the dataset (main list).
    if kind == "main":
        if not sm.is_main:
            return "failed", f"{name} is in the main list but its metadata no longer says every sample is whole blood + bulk. Run `select` again."
        plan: dict[str, list[str] | None] = {"Whole blood": None}
    else:
        plan = dict(sm.bulk_by_type())
        if not plan:
            return "failed", f"{name} has no bulk blood samples in its metadata. Run `select` again."
    if only_types is not None:
        plan = {bt: ids for bt, ids in plan.items() if bt in only_types}
        if not plan:
            return "skipped", f"{name} has no samples of the chosen blood type"

    todo = {bt: ids for bt, ids in plan.items() if force or not config.processed_file(bt, source, project).exists()}
    if not todo:
        return "skipped", f"{name} is already done"

    parts = []
    t0 = time.time()
    try:
        genes = load_genes()
        say("genes", f"{len(genes):,} protein-coding genes ({config.genes_file().name})")
        say("list", f"{name} is in the {kind} list; building {', '.join(config.BLOOD_FOLDER[b] for b in todo)}")

        raw = config.raw_file(source, project)
        got = download(config.raw_url(source, project), raw)
        say("raw", f"{'downloaded' if got else 'already there'}: {raw}  ({raw.stat().st_size / 1e6:.2f} MB)")

        union = None if kind == "main" else sorted({i for ids in todo.values() for i in ids})
        built = build_matrix(raw, genes, union)
        if kind == "main":
            meta_ids = {sid for sid, _, _ in sm.samples}
            if set(built.sample_ids) != meta_ids:
                raise BuildError(f"the raw gene file has {len(built.sample_ids)} samples but the metadata has {len(meta_ids)}; "
                                 f"they must be the same for a main-list dataset")
        say("build", f"{len(genes):,} genes x {len(built.sample_ids):,} samples read")

        n_files = 0
        for bt, ids in todo.items():
            part_built = built if ids is None else built.subset(ids)
            target = config.processed_file(bt, source, project)
            part = target.with_name(target.name + ".part")
            parts.append(part)
            size = write_tsv_gz(part_built.matrix, part)
            info = check_processed(part, genes, part_built.sample_ids, part_built.raw_sample_sums)
            part.replace(target)
            n_files += 1
            say("write", f"{config.BLOOD_FOLDER[bt]:<22} {info['n_samples']:5,} samples  {size / 1e6:6.2f} MB  {target}")
        return "done", f"{name}: {n_files} file(s) in {time.time() - t0:.1f}s"
    except (DownloadError, BuildError, CheckError, ValueError) as e:
        return "failed", f"{name}: {e}"
    finally:
        for part in parts:
            part.unlink(missing_ok=True)      # never leave a half-made file behind


def _types_arg(values: list[str] | None) -> set[str] | None:
    if not values:
        return None
    unknown = [v for v in values if v.lower() not in _TYPE_NAMES]
    if unknown:
        raise SystemExit(f"Unknown blood type {unknown}. Use: {', '.join(sorted(_TYPE_NAMES))}")
    return {_TYPE_NAMES[v.lower()] for v in values}


def _job(args: tuple) -> tuple[str, str, str]:
    source, project, force, types = args
    status, msg = build_one(source, project, force, quiet=True, only_types=types)
    return f"{source.upper()}/{project}", status, msg


def cmd_select(a) -> int:
    main_rows, second_rows = build_lists()
    by = lambda rows: {s: sum(1 for r in rows if r["source"] == s) for s in ("GTEX", "TCGA", "SRA") if any(r["source"] == s for r in rows)}
    print(f"main list   (every sample whole blood + bulk): {len(main_rows)} datasets {by(main_rows)}, "
          f"{sum(r['n_samples'] for r in main_rows):,} samples")
    print(f"second list (every other dataset with a bulk blood sample): {len(second_rows)} datasets {by(second_rows)}, "
          f"{sum(r['n_bulk_blood'] for r in second_rows):,} bulk blood samples")
    print(f"in scope: {len(main_rows) + len(second_rows)} datasets")
    print(f"saved: {config.MAIN_LIST_FILE}\n       {config.SECOND_LIST_FILE}")
    return 0


def cmd_add(a) -> int:
    """Add the datasets found in ANY folder to the lists, leaving the datasets already on the lists untouched."""
    folder = Path(a.source_folder).expanduser()
    if not folder.is_dir():
        print(f"Not a folder: {folder}", file=sys.stderr)
        return 1
    found = find_datasets(folder)
    if not found:
        print(f"No datasets found in {folder}.\nA dataset is recognised by its sample file, for example "
              f"sra.sra.SRP999999.MD.gz, gtex.gtex.BLOOD.MD.gz or tcga.tcga.LAML.MD.gz (in any sub-folder).", file=sys.stderr)
        return 1
    main_rows, second_rows = ([], []) if not config.MAIN_LIST_FILE.exists() else (read_rows(config.MAIN_LIST_FILE), read_rows(config.SECOND_LIST_FILE))
    on_list = {(r["source"].lower(), r["project"]) for r in main_rows + second_rows}

    out = {"main": [], "second": [], "already": [], "no": [], "incomplete": []}
    scratch = Path(tempfile.mkdtemp(prefix="gene_expression_lake_add_")) if a.dry_run else None
    for (source, project), files in found.items():
        name = f"{source.upper()}/{project}"
        if (source, project) in on_list:
            out["already"].append(name)
            continue
        kinds = {f.name.split(".")[1] for f in files}
        needed = ["recount_project"] + (["recount_pred"] if source == "sra" else [])
        missing = [k for k in needed if k not in kinds]
        if missing:
            out["incomplete"].append(f"{name}: missing the {', '.join(missing)} file")
            continue
        if a.dry_run:                      # try it in a temporary place; nothing is copied into the lake
            dest = scratch / source / project
            dest.mkdir(parents=True, exist_ok=True)
            for f in files:
                (dest / f.name).write_bytes(f.read_bytes())
            sm = summarize(source, project, scratch)
        else:
            copy_into_metadata(source, project, files)
            sm = summarize(source, project)
        placed = row_for(sm)
        if placed is None:
            out["no"].append(f"{name}: {refusal_reason(sm)}")
        elif placed[0] == "main":
            out["main"].append(name)
            main_rows.append(placed[1])
        else:
            out["second"].append(f"{name} ({placed[1]['group']})")
            second_rows.append(placed[1])

    def show(title: str, items: list[str]) -> None:
        if items:
            print(f"\n{title}: {len(items)}")
            for x in items[:15]:
                print(f"   {x}")
            if len(items) > 15:
                print(f"   ... {len(items) - 15} more")

    print(f"{len(found)} datasets found in {folder}" + ("  (dry run: nothing was copied or written)" if a.dry_run else ""))
    show("ADDED to the MAIN list (every sample whole blood + bulk)", out["main"])
    show("ADDED to the SECOND list", out["second"])
    show("ALREADY on a list (left as they are)", out["already"])
    show("NOT added: not in scope", out["no"])
    show("NOT added: incomplete metadata", out["incomplete"])
    added = len(out["main"]) + len(out["second"])
    if added and not a.dry_run:
        save_lists(main_rows, second_rows)
        print(f"\nLists saved: {config.MAIN_LIST_FILE.name} ({len(main_rows)} datasets), "
              f"{config.SECOND_LIST_FILE.name} ({len(second_rows)} datasets).")
        print("Next: build them with  build-all  (finished datasets are skipped).")
    elif not added:
        print("\nNothing was added to the lists.")
    return 0


def cmd_list(a) -> int:
    lists = load_lists()
    for which, title in (("main", "MAIN list"), ("second", "SECOND list")):
        keys = [k for k, v in lists.items() if v == which]
        print(f"{title}: {len(keys)} datasets")
        for s, p in keys[: a.top]:
            print(f"   {s.upper():5} {p}")
        if len(keys) > a.top:
            print(f"   ... {len(keys) - a.top} more (use --top N)")
    return 0


def cmd_build(a) -> int:
    status, msg = build_one(a.source, a.project, a.force, only_types=_types_arg(a.blood_type))
    if status == "failed":
        print(f"FAILED {msg}", file=sys.stderr)
        return 1
    _say(status, msg)
    return 0


def cmd_build_all(a) -> int:
    lists = load_lists()
    keys = [k for k, v in lists.items() if a.list in ("all", v)]
    types = _types_arg(a.blood_type)
    if a.limit:
        keys = keys[: a.limit]
    print(f"{len(keys)} datasets selected ({a.list} list).")
    if a.dry_run:
        need = [(s, p) for s, p in keys if not config.raw_file(s, p).exists()]
        total, unknown = 0, 0
        for s, p in need:
            size = remote_size(config.raw_url(s, p))
            if size is None:
                unknown += 1
            else:
                total += size
        print(f"{len(need)} raw gene files still to download: {total / 1e9:.2f} GB"
              + (f"  ({unknown} sizes unknown)" if unknown else ""))
        print("nothing was downloaded or built (dry run).")
        return 0
    load_genes()          # make the gene list once, before the workers start
    log_dir = config.LAKE_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"build_{datetime.now():%Y%m%d_%H%M%S}.log"
    counts = {"done": 0, "skipped": 0, "failed": 0}
    with ProcessPoolExecutor(max_workers=a.workers) as ex, open(log, "w") as lf:
        for i, (name, status, msg) in enumerate(ex.map(_job, [(s, p, a.force, types) for s, p in keys]), 1):
            counts[status] += 1
            line = f"[{i}/{len(keys)}] {status.upper():7} {msg}"
            print(line, flush=True)
            lf.write(line + "\n")
    print(f"\nfinished: {counts['done']} done, {counts['skipped']} skipped, {counts['failed']} failed. Log: {log}")
    return 1 if counts["failed"] else 0


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Gene expression files for Human + bulk RNA-seq + blood datasets.")
    sub = p.add_subparsers(dest="cmd", required=True)
    types_help = "only this blood type (repeatable): " + ", ".join(sorted(_TYPE_NAMES))

    s = sub.add_parser("select", help="make the two lists of datasets in scope, from the metadata files")
    s.set_defaults(fn=cmd_select)

    s = sub.add_parser("add", help="add the datasets found in ANY folder to the lists (existing ones are left alone)")
    s.add_argument("--from", dest="source_folder", required=True, help="a folder that contains the datasets' metadata files, in any layout")
    s.add_argument("--dry-run", action="store_true", help="only show what would happen; nothing is copied or written")
    s.set_defaults(fn=cmd_add)

    s = sub.add_parser("list", help="show the two lists")
    s.add_argument("--top", type=int, default=10)
    s.set_defaults(fn=cmd_list)

    s = sub.add_parser("build", help="build ONE dataset from the lists")
    s.add_argument("--source", required=True, help="sra, gtex or tcga")
    s.add_argument("--project", required=True, help="the dataset id, for example SRP158591")
    s.add_argument("--blood-type", action="append", help=types_help)
    s.add_argument("--force", action="store_true", help="redo it even if it is already done")
    s.set_defaults(fn=cmd_build)

    s = sub.add_parser("build-all", help="build every dataset in the lists")
    s.add_argument("--list", choices=["main", "second", "all"], default="all", help="which list (default all)")
    s.add_argument("--blood-type", action="append", help=types_help)
    s.add_argument("--workers", type=int, default=2, help="datasets built at the same time (default 2)")
    s.add_argument("--limit", type=int, default=None, help="only the first N datasets")
    s.add_argument("--dry-run", action="store_true", help="only show how much would be downloaded")
    s.add_argument("--force", action="store_true", help="redo datasets that are already done")
    s.set_defaults(fn=cmd_build_all)

    a = p.parse_args(argv)
    sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
