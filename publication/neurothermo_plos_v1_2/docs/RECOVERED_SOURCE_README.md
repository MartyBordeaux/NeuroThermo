# Recovered upstream source bundles

The four source trees were reconstructed from the exact Git blobs on `freeze-working-code-through-v1.3.1`. Each bundle is Base85-encoded zlib-compressed JSON mapping relative paths to UTF-8 source text. `RECOVERED_SOURCE_HASHES.tsv` records the immutable Git blob SHA, compressed payload SHA-256, JSON payload SHA-256, and a deterministic restored-tree SHA-256. Historical ZIP SHA-256 values from the freeze manifests are retained as provenance only because the JSON bundle is a different container. No scientific source code was modified during restoration.
