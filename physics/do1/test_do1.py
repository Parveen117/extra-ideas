import unittest
import do1_diagonal_observer as m


class TestDO1(unittest.TestCase):
    def test_all(self):
        out, rec = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
