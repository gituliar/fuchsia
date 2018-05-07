#!/usr/bin/env sage-python
import logging
import os.path
import sys

import fuchsia
fuchsia.USE_MAPLE = True

import sage.all
print "Welcome to the Normalize Assistant!"

def print_usage():
    print "Usage:"
    print "  normalize.py <matrix>"

if len(sys.argv) != 2:
    print "invalid arguments: normalize.py needs exactly 1 argument"
    print_usage()
    exit(1)
fname = sys.argv[1]

if not os.path.isfile(fname):
    print "'%s' does not exist or is not a file" % fname
    exit(1)

fuchsia.logger.setLevel(logging.INFO)
m0 = fuchsia.import_matrix_from_file(fname)
x, ep = sage.all.var("x ep")

print "You are normalizing matrix '%s' in (%s, %s)" % (fname, x, ep)

m = fuchsia.FuchsianSystem.from_M(m0, x, ep)
na = fuchsia.NormalizeAssistant(m)
na.start()
