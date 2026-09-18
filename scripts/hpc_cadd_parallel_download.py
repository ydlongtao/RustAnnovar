#!/usr/bin/env python3
"""Bounded, resumable official CADD downloads; expose only MD5-verified files."""
import argparse
import concurrent.futures
import fcntl
import hashlib
import http.client
import json
import os
import shutil
import time
import urllib.request
from pathlib import Path


def digest(path):
    value = hashlib.md5()
    with path.open("rb") as handle:
        for data in iter(lambda: handle.read(1048576), b""):
            value.update(data)
    return value.hexdigest()


def request(url, start, end, length=None, validator=None):
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}", "Accept-Encoding": "identity"})
    response = urllib.request.urlopen(req, timeout=60)
    try:
        prefix = f"bytes {start}-{end}/"
        if response.status != 206 or not response.headers.get("Content-Range", "").startswith(prefix):
            raise ValueError("invalid CADD HTTP range response")
        actual = int(response.headers["Content-Range"][len(prefix):])
        tag = response.headers.get("ETag") or response.headers.get("Last-Modified")
        if (length is not None and actual != length) or (validator is not None and tag != validator):
            raise ValueError("CADD source identity changed")
        if not tag:
            raise ValueError("CADD source lacks identity validator")
        return response, actual, tag
    except Exception:
        response.close()
        raise


def download(url, destination, checksum, workers, chunk_size):
    if destination.exists():
        if digest(destination) != checksum:
            raise ValueError(f"existing CADD file fails official MD5: {destination}")
        return
    with request(url, 0, 0)[0] as response:
        length = int(response.headers["Content-Range"].split("/")[1])
        validator = response.headers.get("ETag") or response.headers["Last-Modified"]
    if shutil.disk_usage(destination.parent).free < 2 * length + 1073741824:
        raise RuntimeError("insufficient filesystem free space for chunks and merge")
    chunks = destination.with_name(destination.name + ".chunks")
    chunks.mkdir(exist_ok=True)
    manifest = {"url": url, "length": length, "validator": validator, "chunk_size": chunk_size}
    metadata = chunks / "source.json"
    if metadata.exists():
        if json.loads(metadata.read_text()) != manifest:
            raise ValueError("existing chunks belong to another CADD source")
    else:
        metadata.write_text(json.dumps(manifest, indent=2) + "\n")
    old = destination.with_name(destination.name + ".part")

    def worker(index):
        start = index * chunk_size
        end = min(length, start + chunk_size) - 1
        size = end - start + 1
        final = chunks / f"{index:06d}"
        partial = chunks / f"{index:06d}.part"
        if final.exists():
            if final.stat().st_size != size:
                raise ValueError("invalid completed CADD chunk size")
            return
        # Existing single-stream bytes remain intact; reuse complete/partial prefix.
        if not partial.exists() and old.exists() and old.stat().st_size > start:
            available = min(size, old.stat().st_size - start)
            with old.open("rb") as source, partial.open("wb") as target:
                source.seek(start)
                while available:
                    data = source.read(min(1048576, available))
                    if not data:
                        raise ValueError("old CADD prefix was truncated")
                    target.write(data)
                    available -= len(data)
        for attempt in range(10):
            offset = partial.stat().st_size if partial.exists() else 0
            if offset > size:
                raise ValueError("oversized partial CADD chunk")
            if offset == size:
                partial.rename(final)
                print(f"completed chunk {index}", flush=True)
                return
            try:
                response, _, _ = request(url, start + offset, end, length, validator)
                with response, partial.open("ab") as handle:
                    remaining = size - offset
                    while remaining:
                        data = response.read(min(1048576, remaining))
                        if not data:
                            raise OSError("truncated CADD range")
                        handle.write(data)
                        remaining -= len(data)
                    if response.read(1):
                        raise ValueError("oversized CADD range body")
                    handle.flush()
                    os.fsync(handle.fileno())
            except (OSError, TimeoutError, http.client.HTTPException) as error:
                print(f"chunk {index}: attempt {attempt+1}: {error}", flush=True)
                if attempt == 9:
                    raise
                time.sleep(min(30, 2 ** attempt))
        # A successful final transfer still needs atomic completion.
        if partial.stat().st_size != size:
            raise ValueError("incomplete CADD chunk")
        partial.rename(final)

    count = (length + chunk_size - 1) // chunk_size
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(worker, range(count)))
    merged = destination.with_name(destination.name + ".parallel.part")
    with merged.open("wb") as target:
        for index in range(count):
            with (chunks / f"{index:06d}").open("rb") as source:
                shutil.copyfileobj(source, target, 1048576)
        target.flush()
        os.fsync(target.fileno())
    if merged.stat().st_size != length or digest(merged) != checksum:
        raise ValueError("merged CADD file fails official MD5; pieces preserved")
    merged.rename(destination)
    print(f"official MD5 verified: {destination}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--build", choices=["hg19", "hg38"], required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error("workers must be between 1 and 8")
    assembly = "GRCh37" if args.build == "hg19" else "GRCh38"
    directory = args.root / "humandb" / "cadd-v1.7" / assembly
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "download.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        base = f"https://krishna.gs.washington.edu/download/CADD/v1.7/{assembly}"
        with urllib.request.urlopen(f"{base}/MD5SUMs", timeout=60) as response:
            text = response.read(65537).decode()
        if len(text) > 65536:
            raise ValueError("oversized CADD checksums")
        checksums = dict((line.split()[1], line.split()[0]) for line in text.splitlines() if line.strip())
        (directory / "MD5SUMs").write_text(text)
        for name in ["whole_genome_SNVs.tsv.gz.tbi", "whole_genome_SNVs.tsv.gz"]:
            checksum = checksums[name]
            if len(checksum) != 32 or any(c not in "0123456789abcdef" for c in checksum):
                raise ValueError("invalid official CADD checksum")
            download(f"{base}/{name}", directory / name, checksum, args.workers, 64 * 1048576)


if __name__ == "__main__":
    main()
