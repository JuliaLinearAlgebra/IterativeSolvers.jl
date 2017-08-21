# [(Inverse) power method](@id PowerMethod)

Solves the eigenproblem $Ax = λx$ approximately where $A$ is a general linear map. By default converges towards the dominant eigenpair $(λ, x)$ such that $|λ|$ is largest. Shift-and-invert can be applied to target a specific eigenvalue near `shift` in the complex plane.

## Usage

```@docs
powm
powm!
invpowm
invpowm!
```

## Implementation details
Storage requirements are 3 vectors: the approximate eigenvector `x`, the residual vector `r` and a temporary. The residual norm lags behind one iteration, as it is computed when $Ax$ is performed. Therefore the final resdiual norm is even smaller.
