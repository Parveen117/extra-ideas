import unittest
import oa1_once_around_with_the_full_form as m


class TestOA1(unittest.TestCase):
    def test_all(self):
        out, num = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)
        self.assertAlmostEqual(num['per_circuit_mas']/num['GR2_per_circuit_mas'], 1.5, places=3)


if __name__ == '__main__':
    unittest.main()
