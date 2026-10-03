"""Markdown table of a proof replay from its aggregate receipt (replays/receipts/<ID>.jsonl, written by replay_all.py).  The receipt
is validated first with replay_all.read_receipt (canonical records, strict sequence, chained prev_sha256, per-record self_sha256,
schema, a unique final end / abort record); then for every replayed log: family, log sha256 prefix, shards, leaves checked / ok /
failures (from the shard summaries recorded in the receipt), validation problems; the gate results; the final status and digest.
Usage: summarize.py RECEIPT"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from replay_all import read_receipt
def main():
    raw = open(sys.argv[1], 'rb').read()
    recs, digest, problems = read_receipt(raw)
    if problems: print('RECEIPT INVALID:', problems); sys.exit(1)
    start = recs[0]['data']; end = recs[-1]['data'] if recs[-1]['kind'] == 'end' else None
    print(f"run {start['run_id']} mode {start['mode']} shards {start['shards']} pin {start['pin_file']} ({start['pin_count']} files, sha256 {start['pin_sha256'][:12]})")
    print('| log | family | log sha256 | shards ok | leaves checked | ok | failures | validation problems |')
    print('|---|---|---|---|---|---|---|---|')
    for r in recs:
        if r['kind'] != 'replay': continue
        d = r['data']; s = [x for x in d['summaries'] if x]
        ok_sh = sum(1 for c in d['exit_codes'].values() if c == 0)
        print(f"| `{d['log']}` | {d['family']} | {d['log_sha256'][:12]} | {ok_sh}/{d['shards']} | {sum(x['leaves_checked'] for x in s)} | "
              f"{sum(x['ok'] for x in s)} | {sum(x['failures'] for x in s)} | {d['validation_problems']} |")
    for r in recs:
        if r['kind'] == 'gate': print(f"gate {' '.join(r['data']['command'])}: exit {r['data']['exit']}")
    print('end:', {k: end[k] for k in ('status', 'end_state_unchanged', 'failures')} if end else 'NO END RECORD (aborted or incomplete run)')
    print('final digest:', digest)
if __name__ == '__main__':
    main()
