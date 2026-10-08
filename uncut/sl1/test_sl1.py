import importlib.util, os
spec = importlib.util.spec_from_file_location('sl1', os.path.join(os.path.dirname(__file__), 'sl1_unbiased_histories_raise_the_count.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
