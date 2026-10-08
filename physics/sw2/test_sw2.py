import unittest
import sw2_turning_centre_exact as m


class TestSW2(unittest.TestCase):
    def test_all(self):
        out = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
