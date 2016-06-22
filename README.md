*Fuchsia* reduces differential equations for Feynman master integrals to canonical form.

In concrete terms, let's say you have a system of differential equations of this form:

    d/dx F(x, eps) = M(x, eps) F(x, eps)

where `x` is a free variable, `M(x, eps)` is a matrix of rational functions (in `x` and `eps`), and `F(x, eps)` is a column vector of unknown functions we are looking for as a Laurent series in `eps`, an infinitesimal parameter.

The task of *Fuchsia* is to find an equivalent Fuchsian system of this form:

    d/dx G(x, eps) = eps S(x) G(x, eps)

where S(x) = sum_i { S_i / (x - x_i) } and the transformation itself defined by matrix `T(x, eps)' like this:

    F(x, eps) = T(x, eps) G(x, eps)

Such a transformation is useful, because from it one can easily
find J' (and therefore J) as a series in eps.

You can learn about the algorithm *Fuchsia* uses to perform this
transformation from Roman Lee's paper at [1].

*Fuchsia* is available both as a command line utility and as a
(Python) library for SageMath [2]. It will run on most Unix-like
operating systems. You can learn about it's installation and
usage from [3].

[1] https://arxiv.org/abs/1411.0911v1
[2] http://www.sagemath.org/
[3] http://www.gituliar.net/fuchsia/fuchsia.pdf
