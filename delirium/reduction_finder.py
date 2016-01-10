import logging

from   delirium.matrix import coefficient_low, degree_low, dims, new_Matrix, var

log = logging.getLogger('delirium')

class ReductionFinder(object):

    def jordan_form(self, m, x):
        x = var(x)
        mj = m.jordan_form()

    def moser_form(self, m, x):
        x = var(x)
        m0 = coefficient_low(m, x)
        log.info("\nA =\n" + str(m))
        log.info("\nA_0 =\n" + str(m0))
        mm = m
        n = dims(m)[0]
        t = new_Matrix([0]*n**2)
        return mm, t

#        print "\nMatrix valuation:", v, "\n"
#        print "\nMatrix [0] eigenvalues:\n", m0.eigenvalues(), "\n"
#        print "\nMatrix [0] eigenvectors:\n", m0.eigenvectors_left(), "\n"
#        print "\nMatrix [0] Jordan form:\n", m0.jordan_form(transformation=True), "\n"
#        m1 = nb.matrix_coefficient(m, 'x', v+1)
#        l = var('l')
#        print "\nDet [0,1]:\n", (m0/x+m1 - l).det(), "\n"
