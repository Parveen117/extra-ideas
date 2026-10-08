import unittest
import gf1_no_memory_in_every_frame as m


class TestGF1(unittest.TestCase):
    def test_all(self):
        out, forms, turn = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
