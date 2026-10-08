import unittest
import cc1_cost_of_a_content as m


class TestCC1(unittest.TestCase):
    def test_all(self):
        out, num = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
