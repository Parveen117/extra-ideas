import unittest
import rc1_rational_closure as m


class TestRC1(unittest.TestCase):
    def test_all(self):
        out, rec = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
