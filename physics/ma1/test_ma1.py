import unittest
import ma1_mo1_assumption_from_the_law as m


class TestMA1(unittest.TestCase):
    def test_all(self):
        out = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
