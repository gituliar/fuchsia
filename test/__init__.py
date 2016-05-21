import doctest
import time
import unittest

import fuchsia

class TimedTextTestResult(unittest.TextTestResult):
    def startTest(self, test):
        super(unittest.TextTestResult, self).startTest(test)
        test._startTime = time.time()
        self.stream.write(self.getDescription(test))
        self.stream.write(" ... ")
        self.stream.flush()

    def addSuccess(self, test):
        super(unittest.TextTestResult, self).addSuccess(test)
        duration = time.time() - test._startTime
        self.stream.writeln("ok, %.3f sec" % duration)

    def addError(self, test, err):
        super(unittest.TextTestResult, self).addError(test, err)
        duration = time.time() - test._startTime
        self.stream.writeln("ERROR, %.3f sec" % duration)

    def addFailure(self, test, err):
        super(unittest.TextTestResult, self).addFailure(test, err)
        duration = time.time() - test._startTime
        self.stream.writeln("FAIL, %.3f sec" % duration)

    def addSkip(self, test, reason):
        super(unittest.TextTestResult, self).addSkip(test, reason)
        self.stream.writeln("skipped {0!r}".format(reason))

    def addExpectedFailure(self, test, err):
        super(unittest.TextTestResult, self).addExpectedFailure(test, err)
        duration = time.time() - test._startTime
        self.stream.writeln("expected failure, %.3f sec" % duration)

    def addUnexpectedSuccess(self, test):
        super(unittest.TextTestResult, self).addUnexpectedSuccess(test)
        duration = time.time() - test._startTime
        self.stream.writeln("unexpected success, %.3f sec" % duration)

unittest.TextTestRunner.resultclass = TimedTextTestResult

class FastTestLoader(unittest.TestLoader):
    def getTestCaseNames(self, testCaseClass):
        fnNames = super(FastTestLoader, self).getTestCaseNames(testCaseClass)
        return [name for name in fnNames if "slow" not in name]

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(fuchsia))
    return tests

def fast_test_suite():
    loader = FastTestLoader()
    tests = doctest.DocTestSuite(fuchsia)
    tests.addTests(loader.discover("test"))
    return tests

def full_test_suite():
    loader = unittest.TestLoader()
    tests = doctest.DocTestSuite(fuchsia)
    tests.addTests(loader.discover("test"))
    return tests
