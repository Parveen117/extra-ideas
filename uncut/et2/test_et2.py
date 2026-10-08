import importlib.util, os
spec = importlib.util.spec_from_file_location('et2', os.path.join(os.path.dirname(__file__), 'et2_three_partners_zero_surface.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
