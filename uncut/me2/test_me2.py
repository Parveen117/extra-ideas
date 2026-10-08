import importlib.util, os
spec = importlib.util.spec_from_file_location('me2', os.path.join(os.path.dirname(__file__), 'me2_all_planes_and_windows.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
