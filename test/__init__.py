import doctest
import fuchsia

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(fuchsia))
    return tests
