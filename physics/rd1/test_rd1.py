import unittest
import rd1_light_gas_numbers as m


class TestRD1(unittest.TestCase):
    def test_all(self):
        out, num = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
