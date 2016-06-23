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
        sh("bin/fuchsia", "-h")

    def test_fuchsify_henn324(t):
        with Temp() as mfile, Temp() as tfile:
            sh("bin/fuchsia", "fuchsify", "-m", mfile, "-t", tfile,
                    "examples/henn_324.mtx")
            t.assertTransformation("examples/henn_324.mtx", "x", tfile, mfile)
            t.assertIsFuchsian(mfile, "x")

    def test_fuchsify_git409(t):
        with Temp() as mfile, Temp() as tfile:
            sh("bin/fuchsia", "fuchsify", "-m", mfile, "-t", tfile,
                    "examples/git_409.mtx")
            t.assertTransformation("examples/git_409.mtx", "x", tfile, mfile)
            t.assertIsFuchsian(mfile, "x")

    def test_reduce_git410_slow(t):
        with Temp() as mfile, Temp() as tfile:
            sh("bin/fuchsia", "reduce", "--use-maple", "-m", mfile, "-t", tfile,
                    "examples/git_410.mtx")
            t.assertTransformation("examples/git_410.mtx", "x", tfile, mfile)
            t.assertIsReduced(mfile, "x", "eps")

    def test_reduce_lee03_slow(t):
        with Temp() as mfile, Temp() as tfile:
            sh("bin/fuchsia", "reduce", "--use-maple", "-m", mfile, "-t", tfile,
                    "examples/lee_03.mtx")
            t.assertTransformation("examples/lee_03.mtx", "x", tfile, mfile)
            t.assertIsReduced(mfile, "x", "eps")

    def test_reduce_pap01_slow(t):
        with Temp() as mfile, Temp() as tfile:
            sh("bin/fuchsia", "reduce", "--use-maple", "-m", mfile, "-t", tfile,
                    "-f", "m", "-e", "ep", "examples/pap_01.m")
            t.assertTransformation("examples/pap_01.m", "x", tfile, mfile)
            t.assertIsReduced(mfile, "x", "ep")

    def test_reduce_pap02_slow(t):
        with Temp() as mfile, Temp() as tfile:
            sh("bin/fuchsia", "reduce", "--use-maple", "-m", mfile, "-t", tfile,
                    "-f", "m", "-e", "ep", "examples/pap_02.m")
            t.assertTransformation("examples/pap_02.m", "x", tfile, mfile)
            t.assertIsReduced(mfile, "x", "ep")
