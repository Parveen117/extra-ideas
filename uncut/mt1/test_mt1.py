import importlib.util, os
spec = importlib.util.spec_from_file_location('mt1', os.path.join(os.path.dirname(__file__), 'mt1_moving_reader_temperature.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def test_all():
    for k, val in m.out.items():
        if isinstance(val, bool): assert val, k
