import unittest
import or1_one_rule_for_records as m


class TestOR1(unittest.TestCase):
    def test_all(self):
        out, rec = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
