#!/usr/bin/env python3
"""
download_luad_maf.py
====================
Reproducibly download all TCGA-LUAD open-access masked somatic mutation
MAF files from the GDC API and concatenate into a single MAF.gz file.

This reproduces the exact file used in the Sherlock-Lung TCGA analysis:
  - Project: TCGA-LUAD
  - Data type: Masked Somatic Mutation (MuTect2)
  - Access: open
  - Expected: ~618 per-sample MAFs, 558 unique patients, GRCh38

Usage:
    python3 scripts/download_luad_maf.py
    python3 scripts/download_luad_maf.py --output /path/to/output.maf.gz
    python3 scripts/download_luad_maf.py --batch-size 50 --dry-run

Author: Abid Al Reza, PhD
Prepared for: Computational Scientist II, CGR/FNL/NCI (req4512)
"""

import argparse
import gzip
import io
import json
import os
import sys
import tarfile
import time
import urllib.error
import urllib.request

# ---- Configuration ----
GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT  = "https://api.gdc.cancer.gov/data"

GDC_FILTER = {
    "op": "and",
    "content": [
        {
            "op": "=",
            "content": {
                "field": "cases.project.project_id",
                "value": "TCGA-LUAD"
            }
        },
        {
            "op": "=",
            "content": {
                "field": "data_type",
                "value": "Masked Somatic Mutation"
            }
        },
        {
            "op": "=",
            "content": {
                "field": "access",
                "value": "open"
            }
        }
    ]
}

# ---- Helpers ----

def gdc_query_file_ids(verbose=True):
    """Query GDC API for all TCGA-LUAD open-access MAF file IDs."""
    payload = json.dumps({
        "filters": GDC_FILTER,
        "fields":  "file_id,file_name,file_size,cases.submitter_id",
        "format":  "JSON",
        "size":    "1000"
    }).encode("utf-8")

    req = urllib.request.Request(
        GDC_FILES_ENDPOINT,
        data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json"}
    )

    if verbose:
        print("Querying GDC API for TCGA-LUAD MAF files...")

    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
    except urllib.error.URLError as e:
        print(f"ERROR: GDC API query failed: {e}", file=sys.stderr)
        sys.exit(1)

    hits     = data["data"]["hits"]
    file_ids = [h["file_id"] for h in hits]
    total    = data["data"]["pagination"]["total"]

    if verbose:
        print(f"  Found {total} files in GDC ({len(file_ids)} returned)")
        total_size_mb = sum(
            h.get("file_size", 0) for h in hits
        ) / 1_000_000
        print(f"  Total size: {total_size_mb:.1f} MB")

    return file_ids


def gdc_bulk_download_batch(file_ids, batch_num, total_batches, verbose=True):
    """
    POST a list of file IDs to GDC data endpoint.
    Returns raw bytes of the tar.gz response.
    """
    payload = json.dumps({"ids": file_ids}).encode("utf-8")
    req = urllib.request.Request(
        GDC_DATA_ENDPOINT,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept":       "application/json"
        }
    )

    if verbose:
        print(f"  Downloading batch {batch_num}/{total_batches} "
              f"({len(file_ids)} files)...", end="", flush=True)

    retries = 3
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                raw = resp.read()
            if verbose:
                print(f" {len(raw)/1e6:.1f} MB")
            return raw
        except urllib.error.URLError as e:
            if attempt < retries:
                wait = 10 * attempt
                if verbose:
                    print(f"\n  Attempt {attempt} failed ({e}). "
                          f"Retrying in {wait}s...", end="", flush=True)
                time.sleep(wait)
            else:
                print(f"\nERROR: Batch {batch_num} failed after {retries} attempts: {e}",
                      file=sys.stderr)
                raise


def extract_maf_lines_from_tar(raw_bytes):
    """
    Extract all mutation lines from a GDC bulk download tar.gz.
    Returns (header_lines, data_lines) where header_lines is a list
    of comment/column-header lines from the first MAF encountered.
    """
    header_lines = []
    data_lines   = []
    seen_header  = False

    try:
        with tarfile.open(fileobj=io.BytesIO(raw_bytes), mode="r:gz") as tar:
            for member in tar.getmembers():
                if not member.name.endswith(".maf") and \
                   not member.name.endswith(".maf.gz"):
                    continue

                f = tar.extractfile(member)
                if f is None:
                    continue

                raw = f.read()

                # Handle gzipped MAFs inside the tar
                if member.name.endswith(".gz"):
                    raw = gzip.decompress(raw)

                lines = raw.decode("utf-8", errors="replace").splitlines()

                file_header = []
                file_data   = []
                for line in lines:
                    if line.startswith("#") or \
                       (not seen_header and line.startswith("Hugo_Symbol")):
                        file_header.append(line)
                    else:
                        if line.strip():
                            file_data.append(line)

                if not seen_header and file_header:
                    header_lines = file_header
                    seen_header  = True

                data_lines.extend(file_data)

    except tarfile.TarError as e:
        print(f"\n  WARNING: Could not parse tar batch: {e}", file=sys.stderr)

    return header_lines, data_lines


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


# ---- Main ----

def main():
    parser = argparse.ArgumentParser(
        description="Download and concatenate TCGA-LUAD MAFs from GDC API"
    )
    parser.add_argument(
        "--output", "-o",
        default="data/TCGA_LUAD/TCGA_LUAD_somatic.maf.gz",
        help="Output MAF.gz path (default: data/TCGA_LUAD/TCGA_LUAD_somatic.maf.gz)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=100,
        help="Number of files per GDC API batch request (default: 100)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Query GDC and report file count/size without downloading"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Suppress progress output"
    )
    args = parser.parse_args()

    verbose = not args.quiet

    if verbose:
        print("=" * 60)
        print("TCGA-LUAD Somatic MAF Download")
        print("Source: GDC API (open-access masked somatic mutations)")
        print("=" * 60)

    # Step 1: Query file IDs
    file_ids = gdc_query_file_ids(verbose=verbose)

    if not file_ids:
        print("ERROR: No files found. Check GDC API availability.",
              file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        print(f"\nDry run complete. Would download {len(file_ids)} files.")
        print(f"Output would be written to: {args.output}")
        return

    # Step 2: Create output directory
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Step 3: Download in batches and write output
    batches       = list(chunks(file_ids, args.batch_size))
    total_batches = len(batches)
    all_header    = []
    total_lines   = 0
    unique_samples= set()

    if verbose:
        print(f"\nDownloading {len(file_ids)} files in "
              f"{total_batches} batches of {args.batch_size}...")
        print(f"Output: {args.output}\n")

    with gzip.open(args.output, "wt", encoding="utf-8") as out_f:
        for batch_num, batch_ids in enumerate(batches, start=1):
            # Download batch
            raw = gdc_bulk_download_batch(
                batch_ids, batch_num, total_batches, verbose=verbose
            )

            # Extract MAF lines
            header_lines, data_lines = extract_maf_lines_from_tar(raw)

            # Write header only from first batch
            if not all_header and header_lines:
                all_header = header_lines
                for line in header_lines:
                    out_f.write(line + "\n")

            # Write data lines
            for line in data_lines:
                out_f.write(line + "\n")
                # Track unique sample barcodes (column 16 = Tumor_Sample_Barcode)
                parts = line.split("\t")
                if len(parts) > 15:
                    unique_samples.add(parts[15][:12])  # first 12 chars = patient

            total_lines += len(data_lines)

            if verbose:
                print(f"    Cumulative: {total_lines:,} mutations, "
                      f"{len(unique_samples)} unique patients")

            # Brief pause between batches to be polite to GDC API
            if batch_num < total_batches:
                time.sleep(1)

    # Step 4: Summary
    out_size_mb = os.path.getsize(args.output) / 1_000_000

    if verbose:
        print("\n" + "=" * 60)
        print("Download complete.")
        print(f"  Output file:     {args.output}")
        print(f"  File size:       {out_size_mb:.1f} MB")
        print(f"  Total mutations: {total_lines:,}")
        print(f"  Unique patients: {len(unique_samples)}")
        print(f"  Batches:         {total_batches}")
        print("=" * 60)

    if len(unique_samples) < 400:
        print(f"WARNING: Only {len(unique_samples)} unique patients found. "
              f"Expected ~558. Check for download errors.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
