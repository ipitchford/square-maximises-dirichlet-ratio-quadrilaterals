# AI index

Agent-readable map of this archive. Read [STATUS.md](STATUS.md) and [ASSURANCE.md](ASSURANCE.md) before reusing any claim.

## Identity
- Title: *The square maximises the ratio of the first two Dirichlet eigenvalues among quadrilaterals*.
- Version `0.1.0-candidate`, dated 2026-10-03, DOI 10.5281/zenodo.23109314.
- Creator: Anonymous. Publisher: Evidence Press.
- Status: unrefereed candidate. Computer-assisted proof.
- Not formally verified, not independently reproduced, not peer reviewed.

## The statement
- Theorem 1 ([paper/PAPER.md](paper/PAPER.md), §1). For every simple quadrilateral $Q$, convex or not,
  $\lambda_2(Q)/\lambda_1(Q)\le 5/2$ for the Dirichlet Laplacian, with equality only for squares.
- Machine-readable claims, dependencies, external theorems used and exclusions: [CLAIMS.json](CLAIMS.json).

## Where each part of the proof lives
| Claim | Written argument | Certificate data | Verifier / replay | Separately written check |
|---|---|---|---|---|
| Lemmas A, B | [paper/PAPER.md](paper/PAPER.md) §2 | — | — | — |
| Theorem 2 (local step) | §3 (Lemmas 3.1, 3.2, Corollary 3.3) | [local/replay_v2/](local/replay_v2/) | [local/replay_local_v2.sh](local/replay_local_v2.sh) | [verify/xlocal.py](verify/xlocal.py) |
| Theorem 3 (shell) | §4.1–4.4 | [global2/runs/shell0225_cert.jsonl](global2/runs/shell0225_cert.jsonl) | [global2/replay.jl](global2/replay.jl) | [verify/xcover.py](verify/xcover.py) |
| Transport cover | §4.1–4.4 | [global2/runs/transport_cert.jsonl](global2/runs/transport_cert.jsonl) | [global2/replay_t.jl](global2/replay_t.jl) | [verify/xcover.py](verify/xcover.py), [verify/xcert.py](verify/xcert.py) |
| Lemma $B'$ (band) | §4.5 | [bandB/runs/](bandB/runs/) | [bandB/replay.jl](bandB/replay.jl) | [verify/xbandcover.py](verify/xbandcover.py), [verify/xband.py](verify/xband.py) |
| Theorem 5 (darts) | §5 and [paper/SUPPLEMENT_S.pdf](paper/SUPPLEMENT_S.pdf) (source [darts/LEMMAS.md](darts/LEMMAS.md)) | [darts/certs4/](darts/certs4/), [darts/certs/](darts/certs/), [darts/certs2/](darts/certs2/), [darts/manifests/](darts/manifests/), [darts/refcache/](darts/refcache/) | [darts/replay_all.sh](darts/replay_all.sh), pin [darts/review/PINNED.sha256](darts/review/PINNED.sha256) | [verify/xreceipt.py](verify/xreceipt.py), [verify/xdart.py](verify/xdart.py) |

- Soundness specification of every certificate primitive: paper Appendix A.
- Dependency table from theorems to certificates: paper Appendix B.
- Assurance record (reviews, repairs, separately written checks): paper Appendix C.

## Replay
- Entry point: [replay_archive.sh](replay_archive.sh).
- `check` takes minutes: manifests, published-record checks, receipt verification and negative controls.
- The full replays take from 10 minutes to several hours each.
- Results of the fresh-extraction replay of this archive: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md).
- Environment: [ENVIRONMENT.txt](ENVIRONMENT.txt), [requirements.txt](requirements.txt), [environment/julia/](environment/julia/).
- Hashes of every file: [MANIFEST.sha256](MANIFEST.sha256).
- Fresh-replay logs and records: [fresh_replay/](fresh_replay/). Bulk files kept only in the release archive, not in Git:
  [LARGE_FILES.md](LARGE_FILES.md).

## Review and provenance
- Response to the external referee report of 2 October 2026: [REFEREE_RESPONSE.md](REFEREE_RESPONSE.md). The report's text is not
  redistributed.
- Records made during the response: [referee_response/](referee_response/). Examples: the exact Arbon calculation
  [referee_response/arbon_perturbation.py](referee_response/arbon_perturbation.py) and the verified bibliography
  [referee_response/REFERENCES_VERIFIED.md](referee_response/REFERENCES_VERIFIED.md).
- Portability patch of the dart verifier and the fresh proof replay that followed it:
  [darts/review/PORTABILITY_PATCH_2026-10-02.md](darts/review/PORTABILITY_PATCH_2026-10-02.md).
- Prior art: [PRIOR_ART.md](PRIOR_ART.md). Citations: [CITATION_AUDIT.md](CITATION_AUDIT.md), [SOURCES.md](SOURCES.md).
- Provenance (AI-produced; no human mathematical contribution claimed): [PROVENANCE.md](PROVENANCE.md).
- Licences: [LICENSES.md](LICENSES.md).
- Research-gate record: [RESEARCH_GATES.json](RESEARCH_GATES.json).

## Reuse guidance for agents
- Cite Theorem 1 as an unrefereed computer-assisted result, with the DOI.
- Do not describe it as peer reviewed, formally verified or independently reproduced.
- The value $7/3$ for darts is an observation, not a theorem.
- Novelty statements are bounded by [PRIOR_ART.md](PRIOR_ART.md).
