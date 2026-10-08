import unittest
import cs1_count_is_the_source as m


class TestCS1(unittest.TestCase):
    def test_all(self):
        for k, v in m.run().items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
