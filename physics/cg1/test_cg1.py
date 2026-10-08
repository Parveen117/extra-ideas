import unittest
import cg1_cut_graded_generator_on_the_line as m


class TestCG1(unittest.TestCase):
    def test_all(self):
        out, num = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
