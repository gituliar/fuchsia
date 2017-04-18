import os
import subprocess
import tempfile
import unittest

import fuchsia
from   sage.all import SR

class Temp:
    def __enter__(self):
        fd, self.name = tempfile.mkstemp(prefix="fuchsia-")
        os.close(fd)
        return self.name

    def __exit__(self, exc_type, exc_value, traceback):
        os.remove(self.name)

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

    def assertIsReduced(t, m_path, x_name, eps_name):
        M = fuchsia.import_matrix_from_file(m_path)
        x = SR.var(x_name)
        eps = SR.var(eps_name)
        pranks = fuchsia.singularities(M, x).values()
        t.assertEqual(pranks, [0]*len(pranks))
        t.assertTrue(eps not in (M/eps).simplify_rational().variables())

    def test_help(t):
        sh("sage", "-python", "fuchsia.py", "-h")

    def test_fuchsify_henn_324(t):
        with Temp() as mfile, Temp() as tfile:
            sh("sage", "-python", "fuchsia.py", "fuchsify", "-m", mfile, "-t", tfile,
                    "test/data/henn_324.m")
            t.assertTransformation("test/data/henn_324.m", "x", tfile, mfile)
            t.assertIsFuchsian(mfile, "x")

    def test_fuchsify_git_409(t):
        with Temp() as mfile, Temp() as tfile:
            sh("sage", "-python", "fuchsia.py", "fuchsify", "-m", mfile, "-t", tfile,
                    "test/data/git_409.m")
            t.assertTransformation("test/data/git_409.m", "x", tfile, mfile)
            t.assertIsFuchsian(mfile, "x")

    def test_reduce_git_410_slow(t):
        with Temp() as mfile, Temp() as tfile:
            sh("sage", "-python", "fuchsia.py", "reduce", "--use-maple", "-m", mfile, "-t", tfile,
                    "test/data/git_410.m")
            t.assertTransformation("test/data/git_410.m", "x", tfile, mfile)
            t.assertIsReduced(mfile, "x", "eps")

    def test_reduce_lee_3_slow(t):
        with Temp() as mfile, Temp() as tfile:
            sh("sage", "-python", "fuchsia.py", "reduce", "--use-maple", "-m", mfile, "-t", tfile,
                    "test/data/lee_3.m")
            t.assertTransformation("test/data/lee_3.m", "x", tfile, mfile)
            t.assertIsReduced(mfile, "x", "eps")

    def test_reduce_pap_1_slow(t):
        with Temp() as mfile, Temp() as tfile:
            sh("sage", "-python", "fuchsia.py", "reduce", "--use-maple", "-m", mfile, "-t", tfile,
                    "-f", "m", "-e", "ep", "test/data/pap_1.m")
            t.assertTransformation("test/data/pap_1.m", "x", tfile, mfile)
            t.assertIsReduced(mfile, "x", "ep")
