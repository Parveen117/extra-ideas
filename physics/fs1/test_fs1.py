import unittest
import fs1_place_of_the_fine_structure_number as m


class TestFS1(unittest.TestCase):
    def test_all(self):
        out, num = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
