import unittest
import sd1_source_recognized_memory as m


class TestSD1(unittest.TestCase):
    def test_all(self):
        for k, v in m.run().items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
