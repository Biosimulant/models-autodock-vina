"""Input and provenance checks; synthetic command fixtures are not docking evidence."""
import hashlib

import pytest
from biosim import ExecutionContext, ExecutionPolicy
from biosim.signals import make_signal, unwrap_payload
from src.vina_docking_predictor import VinaDockingPredictor

DEFAULTS = {'box_center': [15.19, 53.903, 16.917], 'box_size': [20, 20, 20], 'seed': 2026, 'cpu': 2}


def options(module, value):
    module.set_inputs({'run_options': make_signal(source='test', name='run_options', value=value, emitted_at=0, spec=None)})


def test_partial_overrides_preserve_defaults_and_do_not_leak_previous_overrides():
    model = VinaDockingPredictor(default_run_options=DEFAULTS)
    options(model, {'exhaustiveness': 16, 'seed': 7})
    first = model._resolved_options()
    assert first['box_center'] == DEFAULTS['box_center'] and first['seed'] == 7
    options(model, {'n_poses': 3})
    second = model._resolved_options()
    assert second['seed'] == 2026 and second['n_poses'] == 3 and second['exhaustiveness'] == 8
    assert DEFAULTS['seed'] == 2026


@pytest.mark.parametrize('key,value', [
    ('box_center', [float('nan'), 0, 0]), ('box_center', [float('inf'), 0, 0]),
    ('box_center', [True, 0, 0]), ('box_size', [0, 1, 1]),
    ('box_size', [-1, 1, 1]), ('box_size', [1, float('inf'), 1]),
    ('exhaustiveness', True), ('exhaustiveness', 0), ('n_poses', False),
    ('energy_range', float('nan')), ('energy_range', float('inf')), ('energy_range', True),
    ('cpu', True), ('cpu', -1), ('seed', True), ('seed', 2**31), ('seed', -(2**31)-1),
    ('scoring', 'ad4'),
])
def test_bad_options_rejected_before_launch(key, value):
    model = VinaDockingPredictor(default_run_options=DEFAULTS)
    options(model, {key: value})
    with pytest.raises(ValueError):
        model._resolved_options()


@pytest.mark.parametrize('value', ['invalid', [1, 2]])
def test_non_object_options_rejected(value):
    model = VinaDockingPredictor(default_run_options=DEFAULTS)
    with pytest.raises(ValueError):
        options(model, value)


def test_actual_input_hashes_preserved_on_cli_failure(tmp_path, monkeypatch):
    model = VinaDockingPredictor(default_receptor_pdbqt_path='data/1iep/1iep_receptor.pdbqt', default_ligand_pdbqt_path='data/1iep/1iep_ligand.pdbqt', default_run_options=DEFAULTS, work_dir=str(tmp_path))
    def fail(metadata):
        raise RuntimeError('explicit test failure before CLI launch')
    monkeypatch.setattr(model, '_prepare_runtime', fail)
    result = model.execute({}, context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0, run_end=.01))
    metadata = unwrap_payload(result['run_metadata'])
    assert metadata['status'] == 'error'
    for kind in ['receptor', 'ligand']:
        path = model.model_root / 'data/1iep' / ('1iep_'+kind+'.pdbqt')
        assert metadata['input_sha256'][kind+'_pdbqt'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert metadata['resolved_options']['seed'] == 2026
    assert unwrap_payload(result['pose_summary']) == []


def test_summary_preserves_score_limitations_and_seed():
    model = VinaDockingPredictor(default_run_options=DEFAULTS)
    summary = model._build_docking_summary(pose_records=[{'rank':1,'affinity_kcal_mol':-9.}], options=model._resolved_options(), seed=2026)
    assert summary['seed'] == 2026
    assert 'not measured binding affinities' in summary['score_interpretation']
    assert 'not an experimental structure' in summary['rmsd_reference']

