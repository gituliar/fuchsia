*Fuchsia* reduces differential equations for Feynman master integrals to canonical form.

In concrete terms, let us say we have a system of differential equations of this form:

    ∂F(x,ϵ)/∂x = M(x,ϵ) F(x,ϵ)

where `M(x,ϵ)` is a given matrix of rational functions in `x` and `ϵ`, i.e, a free variable and an infinitesimal parameter.
Our ultimately goal is to find a column vector of unknown functions `F(x,ϵ)` as a Laurent series in `ϵ`, which satisfies our equations.

With the help of *Fuchsia* we can find a transformation matrix `T(x,ϵ)` which turns our system to the equivalent Fuchsian system of this form:

    d/dx G(x,ϵ) = ϵ S(x) G(x,ϵ)

where `S(x) = ∑ᵢ Sᵢ/(x-xᵢ)` and `F(x,ϵ) = T(x,ϵ) G(x,ϵ)`.

Such a transformation is useful, because we can easily solve the equivalent system for `G(x,ϵ)` (see [1]) and then, multiplying it by `T(x,ϵ)`, find `F(x,ϵ)`.

You can learn about the algorithm used in *Fuchsia* to find such transformations from Roman Lee's paper [2].

*Fuchsia* is available both as a command line utility and as a (Python) library for SageMath [3].
It will run on most Unix-like operating systems.

Documentation with more information, installation and usage details is here [4].

  * [1] https://arxiv.org/abs/1304.1806
  * [2] https://arxiv.org/abs/1411.0911
  * [3] http://www.sagemath.org/
  * [4] http://www.gituliar.net/fuchsia/fuchsia.pdf
