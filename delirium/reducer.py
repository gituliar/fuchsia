import logging

from   delirium.matrix import coefficient, degree_low, var

log = logging.getLogger('delirium')

class Reducer(object):

    def reduce_moser(self, m, x):
        x = var(x)
        print "\nInitial matrix A = \n", m, "\n"

        dl = degree_low(m, x)
        m0 = coefficient(m, x, dl)
        print "Leading coeff A_0 =\n", m0, "\n"

#        print "\nMatrix valuation:", v, "\n"
#        print "\nMatrix [0] eigenvalues:\n", m0.eigenvalues(), "\n"
#        print "\nMatrix [0] eigenvectors:\n", m0.eigenvectors_left(), "\n"
#        print "\nMatrix [0] Jordan form:\n", m0.jordan_form(transformation=True), "\n"
#        m1 = nb.matrix_coefficient(m, 'x', v+1)
#        l = var('l')
#        print "\nDet [0,1]:\n", (m0/x+m1 - l).det(), "\n"
