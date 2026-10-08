import unittest
import gm1_charge_on_the_displaced_centre as m


class TestGM1(unittest.TestCase):
    def test_all(self):
        out, num = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)
        self.assertGreater(num['a_over_half_r_s'], 1e40)


if __name__ == '__main__':
    unittest.main()
