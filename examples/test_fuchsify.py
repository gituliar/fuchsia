import unittest

from   fuchsia.matrix import (import_matrix, var)
from   fuchsia.reduction import (fuchsify, is_fuchsian)

eps, x = var("eps, x")

class TestCase(unittest.TestCase):
    def assertFuchsian(test, filename):
        with open('examples/'+filename) as f:
            m = import_matrix(f)
        mm, t = fuchsify(m, x)
        test.assertEqual(is_fuchsian(mm, x), True)

    def test_git_409(test):
        test.assertFuchsian('git_409.mtx')

    def test_git_410(test):
        test.assertFuchsian('git_410.mtx')


if __name__ == "__main__":
    unittest.main(verbosity=2)
