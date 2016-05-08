import unittest

from   sage.all import matrix
import fuchsia
from   fuchsia import (block_triangular_form, import_matrix_from_file)

class Test(unittest.TestCase):
    def test_block_triangular_form_02(t):
        m = matrix([
          [1, 0, 1, 0],
          [0, 1, 0, 1],
          [1, 0, 1, 0],
          [0, 0, 0, 1]
        ])
        mt, tt, b = block_triangular_form(m)
        t.assertEqual(matrix([
          [1, 1, 0, 0],
          [1, 1, 0, 0],
          [0, 0, 1, 0],
          [0, 0, 1, 1]
        ]), mt)
