import json, os, subprocess, sys
def test_all():
    here = os.path.dirname(os.path.abspath(__file__))
    subprocess.check_call([sys.executable, os.path.join(here, 'cr1_chromium_window.py')])
    assert json.load(open(os.path.join(here, 'CR1_RESULT.json')))['checks']['pass']
