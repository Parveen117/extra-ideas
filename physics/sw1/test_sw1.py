import unittest
import sw1_turning_centre_first_order as m


class TestSW1(unittest.TestCase):
    def test_all(self):
        out, num = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)
        self.assertTrue(35 < num['projected_cos_16p84deg'] < 45)


if __name__ == '__main__':
    unittest.main()
