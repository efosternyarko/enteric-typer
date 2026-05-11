#!/usr/bin/env python3
"""Upload assemblies to Pathogenwatch, create a collection, and retrieve
cgMLST + multi-threshold cluster data for E. coli, Salmonella, or Shigella.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import mimetypes
import os
import sys
import time
from pathlib import Path

import requests

DEFAULT_BASE_URL   = "https://next.pathogen.watch"
FASTA_SUFFIXES     = {".fa", ".fasta", ".fna", ".fas", ".fsa"}
TERMINAL_STATUSES  = {"COMPLETE", "FAILEDQC", "FAILED", "ERROR"}


class PathogenwatchError(RuntimeError):
    pass


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pathogenwatch upload + cgMLST + cluster client")
    p.add_argument("--samplesheet",          required=True)
    p.add_argument("--collection-name",      required=True)
    p.add_argument("--thresholds",           default="5,10,20,50",
                   help="Comma-separated cgMLST allele-difference thresholds")
    p.add_argument("--poll-seconds",         type=int, default=60)
    p.add_argument("--max-wait-seconds",     type=int, default=3600)
    p.add_argument("--tree-poll-seconds",    type=int, default=30,
                   help="Seconds between tree-availability polls")
    p.add_argument("--tree-max-wait-seconds",type=int, default=300,
                   help="Max seconds to wait for collection tree (0 = skip tree)")
    p.add_argument("--sample-output",        required=True)
    p.add_argument("--collection-output",    required=True)
    p.add_argument("--summary-output",       required=True)
    p.add_argument("--tree-output",          default=None)
    p.add_argument("--base-url",             default=os.environ.get("PW_API_BASE_URL", DEFAULT_BASE_URL))
    p.add_argument("--api-key",              default=os.environ.get("PW_API_KEY"))
    return p.parse_args()


# ── HTTP helpers ──────────────────────────────────────────────────────────────

def session_for(api_key: str) -> requests.Session:
    s = requests.Session()
    s.headers.update({"X-API-Key": api_key})
    return s


def request_json(session: requests.Session, method: str, url: str, *,
                 params=None, json_body=None, timeout: int = 120,
                 retries: int = 3):
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            r = session.request(method, url, params=params, json=json_body, timeout=timeout)
            r.raise_for_status()
            return r.json() if r.text else None
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as exc:
            last_exc = exc
            if attempt < retries - 1:
                wait = 5 * (2 ** attempt)   # 5 s, 10 s, 20 s
                print(f"  Connection error (attempt {attempt + 1}/{retries}), "
                      f"retrying in {wait}s: {exc}", file=sys.stderr)
                time.sleep(wait)
    raise last_exc  # type: ignore[misc]


def sha1_file(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ── Samplesheet ───────────────────────────────────────────────────────────────

def read_samplesheet(path: Path) -> list[dict]:
    rows = []
    with path.open() as fh:
        for row in csv.DictReader(fh):
            sid   = (row.get("id") or "").strip()
            fasta = Path((row.get("fasta") or "").strip())
            if not sid or not fasta.exists():
                print(f"WARNING: skipping invalid row: {dict(row)}", file=sys.stderr)
                continue
            rows.append({"id": sid, "fasta": fasta})
    if not rows:
        raise PathogenwatchError("No valid rows in samplesheet")
    return rows


# ── Upload ────────────────────────────────────────────────────────────────────

def upload_to_signed_url(upload_url: str, fasta_path: Path) -> None:
    content_type = mimetypes.guess_type(str(fasta_path))[0] or "application/octet-stream"
    payload = gzip.compress(fasta_path.read_bytes())
    r = requests.put(
        upload_url, data=payload,
        headers={"Content-Type": content_type,
                 "Content-Length": str(len(payload)),
                 "Content-Encoding": "gzip"},
        timeout=300,
    )
    r.raise_for_status()


def upload_genome(session: requests.Session, base_url: str,
                  folder_id: int, sample_id: str, fasta_path: Path) -> dict:
    checksum = sha1_file(fasta_path)
    store = request_json(session, "POST", f"{base_url}/api/genomes/store",
                         params={"checksum": checksum, "type": "assembly"})
    if store.get("upload"):
        upload_to_signed_url(str(store["uploadUrl"]), fasta_path)
    genome = request_json(session, "POST", f"{base_url}/api/genomes/create",
                          json_body={"folderId": folder_id, "checksum": checksum,
                                     "name": sample_id})
    return {"sample": sample_id, "fasta": str(fasta_path),
            "checksum": checksum, "id": genome["id"], "uuid": genome["uuid"]}


# ── Polling ───────────────────────────────────────────────────────────────────

def wait_for_processing(session, base_url, uploaded, poll_seconds, max_wait):
    deadline = time.monotonic() + max_wait
    genome_ids = [int(u["id"]) for u in uploaded]
    details_by_uuid: dict[str, dict] = {}
    groups: list[dict] = []

    while time.monotonic() < deadline:
        details_by_uuid = {
            u["uuid"]: request_json(session, "GET", f"{base_url}/api/genomes/details",
                                    params={"id": u["uuid"]})
            for u in uploaded
        }
        groups = request_json(session, "POST", f"{base_url}/api/genomes/group",
                               json_body={"ids": genome_ids}) or []
        if all(str(d.get("status")) in TERMINAL_STATUSES
               for d in details_by_uuid.values()):
            break
        done = sum(1 for d in details_by_uuid.values()
                   if str(d.get("status")) in TERMINAL_STATUSES)
        print(f"  Waiting for processing… ({done}/{len(uploaded)} done)", file=sys.stderr)
        time.sleep(poll_seconds)

    return details_by_uuid, groups


# ── Collection ────────────────────────────────────────────────────────────────

def create_collection(session, base_url, organism_id, genome_ids, name) -> dict:
    return request_json(
        session, "POST", f"{base_url}/api/collections/create",
        json_body={"organismId": organism_id, "genomeIds": genome_ids, "name": name},
    )


def collection_details(session, base_url, collection_uuid) -> dict:
    return request_json(session, "GET", f"{base_url}/api/collections/details",
                        params={"uuid": collection_uuid})


def _trigger_tree_build(session, base_url, collection_uuid) -> None:
    """Best-effort trigger of server-side tree building. Swallows failures."""
    for endpoint, method, body in [
        (f"{base_url}/api/collections/{collection_uuid}/build-tree", "POST", {}),
        (f"{base_url}/api/collections/build-tree",                   "POST", {"uuid": collection_uuid}),
        # Fetching collection details can trigger lazy tree computation on some versions
        (f"{base_url}/api/collections/details?uuid={collection_uuid}", "GET", None),
    ]:
        try:
            if method == "POST":
                r = session.post(endpoint, json=body, timeout=30)
            else:
                r = session.get(endpoint, timeout=30)
            if r.status_code in (200, 201, 202, 204):
                print(f"  Tree build trigger sent ({endpoint})", file=sys.stderr)
                return
        except Exception:
            pass


def wait_for_tree(session, base_url, collection_uuid, tree_output_path: Path | None,
                  poll_seconds: int = 30, max_wait: int = 600) -> bool:
    """Poll until the collection Newick tree is available, triggering build each cycle."""
    if tree_output_path is None:
        return False

    tree_endpoints = [
        f"{base_url}/api/collections/tree?uuid={collection_uuid}",
        f"{base_url}/api/collections/{collection_uuid}/tree",
    ]

    deadline = time.monotonic() + max_wait
    attempt  = 0
    while time.monotonic() < deadline:
        _trigger_tree_build(session, base_url, collection_uuid)
        for endpoint in tree_endpoints:
            try:
                r = session.get(endpoint, timeout=120)
                if r.status_code == 200 and r.text.strip().startswith("("):
                    tree_output_path.write_text(r.text)
                    print(f"  Collection tree written to {tree_output_path} "
                          f"(attempt {attempt + 1})", file=sys.stderr)
                    return True
            except Exception:
                pass
        attempt += 1
        remaining = max(0, int(deadline - time.monotonic()))
        print(f"  Tree not yet ready (attempt {attempt}, "
              f"{remaining}s remaining) — retrying in {poll_seconds}s…", file=sys.stderr)
        time.sleep(min(poll_seconds, remaining + 1))

    print("  Collection tree unavailable after polling — skipping.", file=sys.stderr)
    return False


# ── Clustering ────────────────────────────────────────────────────────────────

def cluster_details_threshold(session, base_url, genome_uuid, threshold: int) -> dict:
    return request_json(session, "GET", f"{base_url}/api/genomes/cluster/details",
                        params={"id": genome_uuid, "threshold": threshold}) or {}


def trigger_recluster(session, base_url, genome_uuid) -> None:
    try:
        request_json(session, "POST", f"{base_url}/api/genomes/cluster/recluster",
                     json_body={"id": genome_uuid})
    except Exception as e:
        print(f"  recluster request failed for {genome_uuid}: {e}", file=sys.stderr)


def cluster_neighbour_labels(cluster_payload: dict, focal_name: str) -> list[str]:
    nodes = cluster_payload.get("nodes") or {}
    if not isinstance(nodes, dict):
        return []
    return sorted({
        str(n.get("label"))
        for n in nodes.values()
        if isinstance(n, dict) and n.get("label") and str(n.get("label")) != focal_name
    })


# ── Output ────────────────────────────────────────────────────────────────────

def build_column_names(thresholds: list[int]) -> list[str]:
    base = [
        "sample", "fasta",
        "pw_status", "pw_species", "pw_organism_id",
        "pw_genome_id", "pw_genome_uuid", "pw_checksum",
        "pw_collection_id", "pw_collection_uuid", "pw_collection_url",
        "pw_cgmlst_st",
        "pw_tree_available", "pw_amr_available",
    ]
    for t in thresholds:
        base += [f"pw_cluster{t}_status", f"pw_cluster{t}_count", f"pw_cluster{t}_labels"]
    return base


def write_sample_tsv(path: Path, rows: list[dict], thresholds: list[int]) -> None:
    columns = build_column_names(thresholds)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> int:
    args    = parse_args()
    api_key = args.api_key
    if not api_key:
        print("ERROR: PW_API_KEY not set. Use --api-key or export PW_API_KEY.", file=sys.stderr)
        return 1

    base_url   = args.base_url.rstrip("/")
    thresholds = [int(t.strip()) for t in args.thresholds.split(",") if t.strip().isdigit()]
    samples    = read_samplesheet(Path(args.samplesheet))
    session    = session_for(api_key)

    print(f"Creating folder '{args.collection_name}'…", file=sys.stderr)
    folder = request_json(session, "POST", f"{base_url}/api/folders/create",
                          json_body={"name": args.collection_name})

    print(f"Uploading {len(samples)} genome(s)…", file=sys.stderr)
    uploaded = []
    for row in samples:
        print(f"  uploading {row['id']}…", file=sys.stderr)
        uploaded.append(upload_genome(session, base_url, int(folder["id"]),
                                      row["id"], row["fasta"]))

    print("Waiting for Pathogenwatch to process genomes…", file=sys.stderr)
    details_by_uuid, groups = wait_for_processing(
        session, base_url, uploaded,
        args.poll_seconds, args.max_wait_seconds,
    )

    supported_groups = [g for g in groups if g.get("supported")]
    top_group = max(supported_groups, key=lambda g: len(g.get("ids", [])), default=None)

    collection:      dict | None = None
    collection_meta: dict        = {}
    tree_downloaded              = False

    if top_group and len(uploaded) > 1:
        organism_id = str(top_group["organismId"])
        genome_ids  = [int(x) for x in top_group["ids"]]
        print(f"Creating collection (organism {organism_id}, {len(genome_ids)} genomes)…",
              file=sys.stderr)
        collection      = create_collection(session, base_url, organism_id,
                                            genome_ids, args.collection_name)
        collection_meta = collection_details(session, base_url, str(collection["uuid"]))

        # Poll until tree is built — triggers build attempts on each cycle
        tree_output = Path(args.tree_output) if args.tree_output else None
        tree_downloaded = wait_for_tree(
            session, base_url, str(collection["uuid"]), tree_output,
            poll_seconds=args.tree_poll_seconds,
            max_wait=args.tree_max_wait_seconds,
        )

    # Per-genome cluster searches at each threshold
    rows: list[dict] = []
    for row in samples:
        uploaded_genome = next(u for u in uploaded if u["sample"] == row["id"])
        genome_uuid     = str(uploaded_genome["uuid"])
        details         = details_by_uuid.get(genome_uuid, {})
        genome_status   = str(details.get("status", "UNKNOWN"))

        result: dict = {
            "sample":              row["id"],
            "fasta":               str(row["fasta"]),
            "pw_status":           genome_status,
            "pw_species":          details.get("species",    "NA"),
            "pw_organism_id":      details.get("organismId", "NA"),
            "pw_genome_id":        uploaded_genome["id"],
            "pw_genome_uuid":      genome_uuid,
            "pw_checksum":         uploaded_genome["checksum"],
            "pw_collection_id":    collection.get("id")   if collection else "NA",
            "pw_collection_uuid":  collection.get("uuid") if collection else "NA",
            "pw_collection_url":   collection.get("url")  if collection else "NA",
            "pw_cgmlst_st":        details.get("cgmlstSt", "NA"),
            "pw_tree_available":   str(tree_downloaded),
            "pw_amr_available":    str(bool(details.get("hasAmr"))),
        }

        if genome_status in {"COMPLETE", "FAILEDQC"}:
            for threshold in thresholds:
                try:
                    cp     = cluster_details_threshold(session, base_url, genome_uuid, threshold)
                    status = str(cp.get("status", "UNKNOWN"))
                    if status != "READY":
                        trigger_recluster(session, base_url, genome_uuid)
                        time.sleep(args.poll_seconds)
                        cp     = cluster_details_threshold(session, base_url, genome_uuid, threshold)
                        status = str(cp.get("status", "UNKNOWN"))
                    labels = cluster_neighbour_labels(cp, str(details.get("name") or row["id"]))
                    result[f"pw_cluster{threshold}_status"] = status
                    result[f"pw_cluster{threshold}_count"]  = len(labels)
                    result[f"pw_cluster{threshold}_labels"] = ";".join(labels) if labels else "NA"
                except Exception as exc:
                    print(f"  Cluster lookup failed for {genome_uuid} "
                          f"at threshold {threshold}: {exc}", file=sys.stderr)
                    result[f"pw_cluster{threshold}_status"] = "ERROR"
                    result[f"pw_cluster{threshold}_count"]  = "NA"
                    result[f"pw_cluster{threshold}_labels"] = "NA"
        else:
            for threshold in thresholds:
                result[f"pw_cluster{threshold}_status"] = "NOT_PROCESSED"
                result[f"pw_cluster{threshold}_count"]  = "NA"
                result[f"pw_cluster{threshold}_labels"] = "NA"

        rows.append(result)

    write_sample_tsv(Path(args.sample_output), rows, thresholds)
    Path(args.collection_output).write_text(
        json.dumps({"folder": folder, "groups": groups,
                    "collection": collection, "collectionMeta": collection_meta}, indent=2) + "\n"
    )
    Path(args.summary_output).write_text(
        json.dumps({"samples": rows, "thresholds": thresholds}, indent=2) + "\n"
    )

    print(f"Done. {len(rows)} samples written to {args.sample_output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except requests.HTTPError as exc:
        print(f"HTTP error: {exc.response.status_code} {exc.response.text[:300]}", file=sys.stderr)
        sys.exit(1)
    except PathogenwatchError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
