import importlib.util, os
spec = importlib.util.spec_from_file_location('qh1', os.path.join(os.path.dirname(__file__), 'qh1_count_step_per_recordable_mode.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
