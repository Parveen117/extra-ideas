import unittest
import th1_three_turns_as_harmonics as m


class TestTH1(unittest.TestCase):
    def test_all(self):
        out, rec = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)


if __name__ == '__main__':
    unittest.main()
