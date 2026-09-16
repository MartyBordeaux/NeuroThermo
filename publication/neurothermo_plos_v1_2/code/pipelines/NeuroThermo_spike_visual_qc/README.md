# NeuroThermo spike visual QC

Standalone, read-only visual audit of spike detection in WT and SCA3 current-clamp ABFs.
It does not import, edit, or overwrite the phenotype pipeline.

## Run

```bash
source /root/venv/bin/activate
python3 -m pip install -r requirements.txt \
  --proxy socks5h://127.0.0.1:1080 --timeout 120 --retries 10

nohup env \
  WT_ROOT="/root/neurothermo/WT" \
  SCA3_ROOT="/root/neurothermo/SCA3" \
  OUTPUT="results_spike_visual_qc" \
  bash run_visual_qc.sh > spike_visual_qc.log 2>&1 &

echo $! > spike_visual_qc.pid
tail -f spike_visual_qc.log
```

## Reading the output

- `02_trace_pages/`: every sweep, with primary and threshold-sensitive events marked.
- `03_event_zooms/ambiguous_event_zooms.pdf`: enlarged waveforms around a balanced subset of disputed events.
- `01_tables/spike_candidates.csv`: every candidate peak and its decision under all detector combinations.
- `01_tables/detector_counts.csv`: counts and rates for every sweep and detector.
- `04_fi_visual/`: raw cell F-I curves, detector comparison, and peak-height distributions. No p-values are used.

Markers on trace pages:

- red filled circle: primary detector, prominence 10 mV and height -20 mV;
- orange open circle: primary event rejected only by the strict -10 mV height;
- blue x: accepted by relaxed 8/-30 but rejected by the primary detector;
- black dot: primary event also retained by 10/-10.

The horizontal dotted lines are -30, -20, and -10 mV. Inspect orange events in WT first: they are the events responsible for much of the detector-dependent change.
