import unittest
import subprocess
import tempfile

import fuchsia
from   sage.all import SR

def sh(*cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        raise Exception("Command %s exited with code %s" % (cmd, p.returncode))
    return stdout, stderr

class Test(unittest.TestCase):
    def assertTransformation(t, m1_path, x_name, t_path, m2_path):
        M1 = fuchsia.import_matrix_from_file(m1_path)
        T = fuchsia.import_matrix_from_file(t_path)
        M2 = fuchsia.import_matrix_from_file(m2_path)
        t.assertEqual(M2.simplify_rational(),
                fuchsia.transform(M1, SR.var(x_name), T).simplify_rational())

    def assertIsFuchsian(t, m_path, x_name):
        M = fuchsia.import_matrix_from_file(m_path)
        x = SR.var(x_name)
        pranks = fuchsia.singularities(M, x).values()
        t.assertEqual(pranks, [0]*len(pranks))

    def test_help(t):
        sh("bin/fuchsia", "-h")

    def test_fuchsify_1(t):
        with tempfile.NamedTemporaryFile() as mf:
            with tempfile.NamedTemporaryFile() as tf:
                sh("bin/fuchsia", "fuchsify", "-m", mf.name, "-t", tf.name,
                        "examples/git_409.mtx")
                t.assertTransformation("examples/git_409.mtx", "x", tf.name, mf.name)
                t.assertIsFuchsian(mf.name, "x")
