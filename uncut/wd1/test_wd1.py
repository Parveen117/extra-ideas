import importlib.util, os
spec = importlib.util.spec_from_file_location('wd1', os.path.join(os.path.dirname(__file__), 'wd1_whole_turns_only_in_three.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
