# Provenance

- **Creator:** Anonymous (Evidence Press standing protocol). **Publisher:** Evidence Press.
- **Human role:**
  - On 30 September 2026 the publisher chose the problem ("problem 12") from a list of open problems in mathematical physics compiled
    by AI agents, and asked for a proof.
  - Later instructions concerned scheduling, compute limits and scope. Examples: "Use research-run as the launcher"; "Don't run an
    independent dart check".
  - On 2 October 2026 the publisher supplied an external referee report and asked for it to be actioned in full and for the
    candidate to be published.
  - No human mathematical contribution is claimed.
- **AI roles:**
  - A research agent (Anthropic Claude) designed the proof, wrote all code, ran the certified covers and replays, and wrote the paper.
  - Separate Claude agents carried out these tasks:
    - wrote the separately written checkers in `verify/`;
    - verified the bibliography (`referee_response/REFERENCES_VERIFIED.md`);
    - typeset Supplement S with a token-level identity check (`paper/supplement_check.py`);
    - carried out the bounded prior-art searches.
  - An OpenAI model ("Sol"), accessed through the API, reviewed the arguments and code adversarially in several rounds. The
    findings and repairs are summarised in paper Appendix C. Review transcripts are not redistributed.
  - The external referee report of 2 October 2026 was supplied by the publisher. Its author is not recorded. Its text and audit
    bundle are not redistributed; each item and its disposition are in [REFEREE_RESPONSE.md](REFEREE_RESPONSE.md).
- **Compute:**
  - One Apple M2 machine (8 cores, 16 GB, macOS 27.0), shared with other workloads.
  - Long jobs ran through a local admission queue that reserves CPU and memory.
  - The largest jobs were the dart recertification and proof replays (several hours each) and the transport cover and its replay
    (about 2 h each on 3 threads).
- **Third-party material:** none redistributed. Cited works are referenced, not included.
