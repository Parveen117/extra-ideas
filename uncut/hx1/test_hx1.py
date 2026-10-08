import importlib.util, os
spec = importlib.util.spec_from_file_location('hx1', os.path.join(os.path.dirname(__file__), 'hx1_fall_is_a_turn_of_cuts.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
