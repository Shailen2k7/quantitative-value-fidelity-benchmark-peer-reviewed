# Quantitative Value Fidelity Benchmark (pilot)

Open harness, corpus, and frozen outputs for the paper "Do large language models
preserve the numbers? A benchmark, an open evaluation harness, and a single model
pilot of quantitative value fidelity in patient facing medical simplification."

## Contents
- corpus.json        Synthetic clinician-style source statements with gold-annotated values and severity weights
- outputs.json       Frozen model simplifications under two prompt conditions (plain, preserve-numbers)
- harness.py         Deterministic extraction, normalisation, alignment, classification, and metrics (Wilson CIs)
- results.json       Computed metrics for both conditions
- audit_plain.csv    Per-value classification audit, plain condition
- audit_preserve.csv Per-value classification audit, preserve-numbers condition

## Reproduce
    python3 harness.py
Outputs results.json and the two audit CSVs.

## License
Released for research reproducibility. Add a license of your choice (e.g. MIT for code, CC-BY for data) before publishing.
