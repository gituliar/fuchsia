*Fuchsia* reduces differential equations for Feynman master integrals to canonical form.

In concrete terms, let us say we have a system of differential equations of this form:

    ∂f(x,ϵ)/∂x = 𝕄(x,ϵ) f(x,ϵ)

where `𝕄(x,ϵ)` is a given matrix of rational functions in `x` and `ϵ`, i.e, a free variable and an infinitesimal parameter.
Our ultimately goal is to find a column vector of unknown functions `f(x,ϵ)` as a Laurent series in `ϵ`, which satisfies our equations.

With the help of *Fuchsia* we can find a transformation matrix `𝕋(x,ϵ)` which turns our system to the equivalent Fuchsian system of this form:

    ∂g(x,ϵ)/∂x = ϵ 𝕊(x) g(x,ϵ)

where `𝕊(x) = ∑ᵢ 𝕊ᵢ/(x-xᵢ)` and `f(x,ϵ) = 𝕋(x,ϵ) g(x,ϵ)`.

Such a transformation is useful, because we can easily solve the equivalent system for `g(x,ϵ)` (see [1]) and then, multiplying it by `𝕋(x,ϵ)`, find `f(x,ϵ)`.

You can learn about the algorithm used in *Fuchsia* to find such transformations from Roman Lee's paper [2].

*Fuchsia* is available both as a command line utility and as a (Python) library for SageMath [3].
It will run on most Unix-like operating systems.

Documentation with more information, installation and usage details is here [4].

  * [1] https://arxiv.org/abs/1304.1806
  * [2] https://arxiv.org/abs/1411.0911
  * [3] http://www.sagemath.org/
  * [4] http://www.gituliar.net/fuchsia/fuchsia.pdf
