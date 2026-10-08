import unittest
import cd1_closure_and_dimension as m


class TestCD1(unittest.TestCase):
    def test_all(self):
        for k, v in m.run().items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
