import importlib.util, os
spec = importlib.util.spec_from_file_location('xu2', os.path.join(os.path.dirname(__file__), 'xu2_expansion_plane_and_sign_test.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
