import importlib.util, os
spec = importlib.util.spec_from_file_location('ts2', os.path.join(os.path.dirname(__file__), 'ts2_share_of_a_turning_plane.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
