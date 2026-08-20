# Loader for the exact v1.0.1 analysis source, split into readable fragments solely to keep repository writes manageable.
from pathlib import Path
_parts = sorted(Path(__file__).with_name('_analysis_parts').glob('part*.pyfrag'))
_src = ''.join(p.read_text(encoding='utf-8') for p in _parts)
exec(compile(_src, str(Path(__file__).with_name('analysis_original.py')), 'exec'), globals(), globals())
