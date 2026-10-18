"""Verify source manifest, and optionally the separately restored data arrays."""
from pathlib import Path
import argparse
import csv
import hashlib
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', action='store_true', help='also verify restored archived data')
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
manifests = ['PACKAGE_MANIFEST.csv'] + (['DATA_MANIFEST.csv'] if args.data else [])
checked = 0
errors = []
for name in manifests:
    with (root/'docs'/name).open(encoding='utf-8-sig', newline='') as stream:
        for row in csv.DictReader(stream):
            target = root/row['path']
            if not target.is_file():
                errors.append(f'Missing: {row["path"]}')
                continue
            digest = hashlib.sha256()
            with target.open('rb') as data:
                for block in iter(lambda: data.read(1024*1024), b''):
                    digest.update(block)
            if target.stat().st_size != int(row['bytes']) or digest.hexdigest() != row['sha256']:
                errors.append(f'Changed: {row["path"]}')
            checked += 1
if errors:
    print('\n'.join(errors),file=sys.stderr)
    sys.exit(1)
print(f'Verified {checked} files against SHA-256 manifests.')
