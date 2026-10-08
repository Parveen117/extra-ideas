import unittest
import se1_stretch_law as m


class TestSE1(unittest.TestCase):
    def test_all(self):
        for k, v in m.run().items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
