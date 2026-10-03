"""Write a pin list ('sha256  path' lines, sorted) of everything a proof replay depends on: every transitive decision source of
every stage of replay_all.py (bootstrap, replays, merge, gates), the launcher replay_all.sh, the lint allowlist, every manifest,
the old logs (certs/, certs2/), the chart caches (refcache/*: unpickled by the verifier, so pinned and checked before loading)
and the certs3 logs read by the containment tests.  Usage: make_pin.py OUT"""
import sys, os, glob
HERE = os.path.dirname(os.path.realpath(__file__)); sys.path.insert(0, HERE); os.chdir(HERE)
from binding import closure, sha_file
import importlib.util
spec = importlib.util.spec_from_file_location('replay_all_mod', os.path.join(HERE, 'replay_all.py'))
ra = importlib.util.module_from_spec(spec); spec.loader.exec_module(ra)
files = (set(closure(ra.ENTRY_POINTS)) | {'replay_all.sh', 'rigor_allowlist.json'} | set(glob.glob('manifests/*.json'))
         | set(glob.glob('certs/*.jsonl')) | set(glob.glob('certs2/*.jsonl')) | set(glob.glob('refcache/*'))
         | {'certs3/pp.merged.jsonl', 'certs3/sp_spk.merged.jsonl', 'certs3/fp.merged.jsonl'})
out = sys.argv[1]
with open(out, 'w') as fo:
    for f in sorted(files): fo.write(f'{sha_file(f)}  {f}\n')
print(len(files), 'files pinned ->', out)
