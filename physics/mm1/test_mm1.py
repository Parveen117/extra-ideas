import unittest
import mm1_one_memory as m


class TestMM1(unittest.TestCase):
    def test_all(self):
        for k, v in m.run().items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
