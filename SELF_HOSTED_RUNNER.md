# Self-hosted runner setup

The runner executes jobs on the VPS. The frozen benchmark remains on the VPS
and is never committed to Git or uploaded as a workflow artifact.

## One-time server setup

On GitHub, open `MartyBordeaux/NeuroThermo` → **Settings** → **Actions** →
**Runners** → **New self-hosted runner**, choose Linux x64, and copy the
current registration command. GitHub issues a short-lived registration token;
do not store it in source code or send it in chat.

On the VPS, run the following as the normal account that owns the data, replacing
`<registration-command-from-github>` with GitHub's command:

```bash
mkdir -p "$HOME/actions-runner" "$HOME/neurothermo_data"
cd "$HOME/actions-runner"
curl -o actions-runner.tar.gz -L \
  https://github.com/actions/runner/releases/latest/download/actions-runner-linux-x64-2.330.0.tar.gz
tar xzf actions-runner.tar.gz
./config.sh --url https://github.com/MartyBordeaux/NeuroThermo \
  --token <token-from-github> --labels neurothermo-vps --unattended
sudo ./svc.sh install
sudo ./svc.sh start
sudo ./svc.sh status
```

GitHub's page is authoritative for the current runner download URL and exact
registration command. If it differs from the example, use GitHub's command.

## Protected benchmark

Install the benchmark once, with read permission limited to the runner account:

```bash
install -D -m 600 \
  /path/to/frozen_v2_w20_observations.csv \
  "$HOME/neurothermo_data/frozen_v2_w20_observations.csv"
```

The workflow uses that path by default. If a different protected location is
needed, define the repository variable `NEUROTHERMO_BENCHMARK` to its absolute
path in GitHub **Settings** → **Secrets and variables** → **Actions** →
**Variables**.

## Launching a run

Open **Actions** → **Per-cell pipeline** → **Run workflow**. Select `smoke`
first. The runner must show `Idle` with label `neurothermo-vps` before dispatch.
After a successful smoke run, select `preliminary`. Use `resume=true` only for
an interrupted run with the same immutable configuration and output directory.

The workflow checks out a fresh worktree for every job. Results are retained as
a private GitHub Actions artifact for 30 days; raw benchmark inputs are not
uploaded.
