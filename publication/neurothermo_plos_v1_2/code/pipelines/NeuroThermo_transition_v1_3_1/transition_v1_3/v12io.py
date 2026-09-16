from __future__ import annotations
from pathlib import Path
import io, json, zipfile
import pandas as pd

class V12Results:
    def __init__(self, directory=None, archive=None):
        d = Path(directory).expanduser().resolve() if directory else None
        a = Path(archive).expanduser().resolve() if archive else None
        self.directory = d if d and d.exists() else None
        self.archive = a if a and a.exists() else None
        if self.directory is None and self.archive is None:
            raise FileNotFoundError(f"Neither v1.2 result directory nor archive exists: directory={d}, archive={a}")
        self._zip = zipfile.ZipFile(self.archive) if self.archive is not None and self.directory is None else None
        if self._zip is not None:
            roots = {n.split('/')[0] for n in self._zip.namelist() if '/' in n}
            if 'results_transition_ensemble_v1_2' in roots:
                self.root = 'results_transition_ensemble_v1_2'
            elif len(roots) == 1:
                self.root = next(iter(roots))
            else:
                raise ValueError('Cannot identify root directory inside v1.2 archive')
        else:
            self.root = None

    @property
    def source_description(self):
        return str(self.directory if self.directory is not None else self.archive)

    def _zip_name(self, rel):
        return f"{self.root}/{rel}"

    def exists(self, rel):
        if self.directory is not None:
            return (self.directory / rel).exists()
        return self._zip_name(rel) in set(self._zip.namelist())

    def read_csv(self, rel, **kwargs):
        if self.directory is not None:
            return pd.read_csv(self.directory / rel, **kwargs)
        with self._zip.open(self._zip_name(rel)) as f:
            return pd.read_csv(f, **kwargs)

    def read_json(self, rel):
        if self.directory is not None:
            return json.loads((self.directory / rel).read_text(encoding='utf-8'))
        with self._zip.open(self._zip_name(rel)) as f:
            return json.loads(f.read().decode('utf-8'))

    def checkpoint_ids(self):
        if self.directory is not None:
            out = []
            for p in (self.directory/'checkpoints').glob('scenario_*.csv.gz'):
                out.append(int(p.name.split('_')[1].split('.')[0]))
            return sorted(out)
        prefix = self._zip_name('checkpoints/scenario_')
        out=[]
        for n in self._zip.namelist():
            if n.startswith(prefix) and n.endswith('.csv.gz'):
                out.append(int(n.rsplit('/',1)[-1].split('_')[1].split('.')[0]))
        return sorted(out)

    def read_checkpoint(self, scenario_id):
        rel = f"checkpoints/scenario_{int(scenario_id):04d}.csv.gz"
        if self.directory is not None:
            return pd.read_csv(self.directory / rel, compression='gzip')
        with self._zip.open(self._zip_name(rel)) as f:
            raw = f.read()
        return pd.read_csv(io.BytesIO(raw), compression='gzip')

    def close(self):
        if self._zip is not None:
            self._zip.close()
