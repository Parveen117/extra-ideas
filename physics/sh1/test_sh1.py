import unittest
import sh1_capacity_of_a_shell as m


class TestSH1(unittest.TestCase):
    def test_all(self):
        out, rec = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
