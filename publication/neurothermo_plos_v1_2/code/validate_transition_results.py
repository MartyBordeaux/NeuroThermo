#!/usr/bin/env python3
from pathlib import Path
import hashlib, json, math, sys

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'

EXPECTED_ARCHIVE_SHA = '621692879a4296507e662f84110c80917bf09482b3355e7e56768a86f1ff2320'
EXPECTED = {
    'wt_exit': 0.1358293233470019,
    'balance': 0.5,
    'sca3_entry': 0.7978563373093712,
    'active_wt_exit': 0.1483293470542087,
    'active_sca3_entry': 0.6832817402871161,
}

def fail(msg):
    print('FAIL ' + msg, file=sys.stderr)
    raise SystemExit(1)

def ok(msg):
    print('PASS ' + msg)

def load(stage):
    p = DATA / stage / 'RUN_SUMMARY.json'
    if not p.exists(): fail(f'missing {p.relative_to(ROOT)}')
    return json.loads(p.read_text(encoding='utf-8'))

def assert_eq(actual, expected, label):
    if actual != expected: fail(f'{label}: {actual!r} != {expected!r}')
    ok(f'{label}={actual!r}')

def assert_close(actual, expected, label, tol=1e-12):
    if not math.isclose(float(actual), float(expected), rel_tol=0.0, abs_tol=tol):
        fail(f'{label}: {actual!r} != {expected!r}')
    ok(f'{label}={float(actual):.15g}')

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20), b''): h.update(chunk)
    return h.hexdigest()

def check_recorded_hash(summary, key, path, stage):
    exp=summary.get('frozen_sha256',{}).get(key)
    if not exp: fail(f'{stage} does not record frozen SHA for {key}')
    if not path.exists(): fail(f'missing comparable frozen source {path.relative_to(ROOT)}')
    got=sha256(path)
    if got != exp: fail(f'{stage} frozen SHA mismatch {key}: {got} != {exp}')
    ok(f'{stage} frozen SHA {key}={got}')

prov = DATA/'upstream_server_bundle'/'ENDPOINT_TRANSITION_ARCHIVE_SHA256.txt'
if not prov.exists(): fail('missing endpoint-transition archive provenance')
prov_sha=prov.read_text(encoding='utf-8').split()[0]
assert_eq(prov_sha, EXPECTED_ARCHIVE_SHA, 'endpoint-transition archive SHA256')

s0=load('transition_v1_0_results')
s11=load('transition_v1_1_results')
s12=load('transition_v1_2_results')
s121=load('transition_v1_2_1_results')
s13=load('transition_v1_3_results')

assert_eq(s0['version'],'1.0.0','v1.0 version')
assert_eq(s0['biological_pairs'],72,'v1.0 biological pairs')
assert_eq(s0['n_scenarios'],988,'v1.0 scenarios')
assert_eq(s0['n_grid'],41,'v1.0 path grid')
assert_eq(s0['n_path_families'],3,'v1.0 path families')
assert_eq(s0['n_path_rows'],988*3*41,'v1.0 path rows')
cp0=DATA/'transition_v1_0_results'/'checkpoints'
if cp0.is_dir(): assert_eq(len(list(cp0.glob('*.json.gz'))),988*3,'v1.0 checkpoint files')

assert_eq(s11['version'],'1.1.0','v1.1 version')
assert_eq(s11['biological_pairs'],72,'v1.1 biological pairs')
assert_eq(s11['support_scenarios'],988,'v1.1 support scenarios')
assert_eq(s11['new_HR_simulations'],0,'v1.1 new HR simulations')
assert_eq(bool(s11['audit_pass']),True,'v1.1 audit pass')
assert_eq(s11['path_rows'],s0['n_path_rows'],'v1.1 inherited path rows')
assert_close(s11['isi_wt_exit_A'],EXPECTED['wt_exit'],'v1.1 ISI WT-exit')
assert_close(s11['isi_sca3_entry_A'],EXPECTED['sca3_entry'],'v1.1 ISI SCA3-entry')
assert_close(s11['active_wt_exit_A'],EXPECTED['active_wt_exit'],'v1.1 active WT-exit')
assert_close(s11['active_sca3_entry_A'],EXPECTED['active_sca3_entry'],'v1.1 active SCA3-entry')

assert_eq(s12['version'],'1.2.0','v1.2 version')
assert_eq(s12['n_biological_pairs'],72,'v1.2 biological pairs')
assert_eq(s12['core_secure_pairs'],32,'v1.2 core-secure pairs')
assert_eq(s12['n_scenarios'],988,'v1.2 scenarios')
assert_eq(s12['n_intrinsic'],31,'v1.2 intrinsic grid')
assert_eq(s12['n_drive'],31,'v1.2 drive grid')
assert_eq(s12['n_state_rows'],988*31*31,'v1.2 state rows')
cp12=DATA/'transition_v1_2_results'/'checkpoints'
if cp12.is_dir(): assert_eq(len(list(cp12.glob('scenario_*.csv.gz'))),988,'v1.2 scenario checkpoint files')
for k,jk in [('WT_exit','wt_exit'),('balance','balance'),('SCA3_entry','sca3_entry')]:
    assert_close(s12['primary_ISI_boundaries'][k],EXPECTED[jk],f'v1.2 ISI {k}')

assert_eq(s121['version'],'1.2.1','v1.2.1 version')
assert_eq(s121['n_biological_pairs'],72,'v1.2.1 biological pairs')
assert_eq(s121['core_secure_pairs'],32,'v1.2.1 core-secure pairs')
assert_eq(s121['n_scenarios'],988,'v1.2.1 scenarios')
assert_eq(s121['new_HR_simulations'],0,'v1.2.1 new HR simulations')
assert_eq(s121['source_v1_2_state_rows'],s12['n_state_rows'],'v1.2.1 source state rows')
for k,jk in [('wt_exit','wt_exit'),('balance','balance'),('sca3_entry','sca3_entry')]:
    assert_close(s121['primary_ISI_boundaries'][k],EXPECTED[jk],f'v1.2.1 ISI {k}')

assert_eq(s13['version'],'1.3.0','v1.3 numerical version')
assert_eq(s13['combined_surface_reused_from'],'1.2.0','v1.3 reused combined surface')
assert_eq(s13['n_biological_pairs'],72,'v1.3 biological pairs')
assert_eq(s13['core_secure_pairs'],32,'v1.3 core-secure pairs')
assert_eq(s13['n_selected_scenarios'],988,'v1.3 scenarios')
assert_eq(s13['n_intrinsic'],31,'v1.3 intrinsic grid')
assert_eq(s13['n_component'],31,'v1.3 component grid')
assert_eq(sorted(s13['new_surface_modes']),sorted(['kappa_only','J_only']),'v1.3 new surface modes')
assert_eq(s13['new_state_rows'],988*31*31*2,'v1.3 new state rows')
for k,jk in [('wt_exit','wt_exit'),('balance','balance'),('sca3_entry','sca3_entry')]:
    assert_close(s13['primary_ISI_boundaries'][k],EXPECTED[jk],f'v1.3 ISI {k}')

# Cross-check the frozen hashes embedded by the historical runs against the
# exact endpoint/v1.1 files now stored in the publication release.
endpoint = DATA/'endpoint_ensemble_v1_0_results'
v11 = DATA/'transition_v1_1_results'
for key, path in [
    ('endpoint_cells_full_observable.csv', endpoint/'endpoint_cells_full_observable.csv'),
    ('transition_ready_endpoint_support.csv', endpoint/'transition_ready_endpoint_support.csv'),
    ('transition_core_transform.csv', endpoint/'transition_core_transform.csv'),
    ('endpoint_geometry.csv', endpoint/'endpoint_geometry.csv'),
]:
    check_recorded_hash(s0,key,path,'v1.0')

for key, path in [
    ('endpoint_cells_full_observable.csv', endpoint/'endpoint_cells_full_observable.csv'),
    ('transition_ready_endpoint_support.csv', endpoint/'transition_ready_endpoint_support.csv'),
    ('PRIMARY_ISI_STAGING.csv', v11/'PRIMARY_ISI_STAGING.csv'),
    ('staging_boundary_definitions_v1_1.csv', v11/'staging_boundary_definitions_v1_1.csv'),
    ('transition_projection_reference_v1_1.csv', v11/'transition_projection_reference_v1_1.csv'),
    ('transition_projection_transform_v1_1.csv', v11/'transition_projection_transform_v1_1.csv'),
]:
    check_recorded_hash(s12,key,path,'v1.2')

# One frozen v1.0/v1.2 input still requires exact recovery from the server.
q75_expected=s12.get('frozen_sha256',{}).get('cell_q75_protocol_anchors.csv')
assert_eq(q75_expected,'fe986b58f072ef0f00d252ca799e1f96b6e5b998a9acc3eade6a78a7967ab136','recorded q75-anchor SHA256')
q75_candidates=[
    DATA/'transition_v1_0_frozen'/'cell_q75_protocol_anchors.csv',
    DATA/'transition_v1_2_frozen'/'cell_q75_protocol_anchors.csv',
]
found=[p for p in q75_candidates if p.exists()]
if found:
    for p in found:
        if sha256(p)!=q75_expected: fail(f'q75-anchor SHA mismatch in {p.relative_to(ROOT)}')
        ok(f'exact q75-anchor recovered: {p.relative_to(ROOT)}')
else:
    print('PENDING exact cell_q75_protocol_anchors.csv not yet recovered; expected SHA256='+q75_expected)

print('TRANSITION_RESULT_INTEGRITY_PASS')
