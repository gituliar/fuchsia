import unittest

from   delirium.expression import alphabet, new_Expression, var

x = var("x")

class Test(unittest.TestCase):
    def test_alphabet(t):
        ex = new_Expression("1/x^2 + 1/(1+x)")
        t.assertEqual(
            alphabet(ex, x),
            set([x, x+1]))

if __name__ == "__main__":
    unittest.main(verbosity=2)
