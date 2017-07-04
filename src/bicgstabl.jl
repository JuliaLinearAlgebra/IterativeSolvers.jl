export bicgstabl, bicgstabl!

bicgstabl(A, b, l::Int = 2; kwargs...) = bicgstabl!(zeros(b), A, b, l; initial_zero = true, kwargs...)

function bicgstabl!(x, A, b, l::Int = 2;
    Pl = Identity(),
    max_mv_products = min(30, size(A, 1)),
    initial_zero = false,
    tol = sqrt(eps(real(eltype(b)))),
    visitor = (a...) -> nothing
)
    T = eltype(b)
    n = size(A, 1)

    mv_products = 0
    residuals = real(T)[]

    rs = zeros(T, n, l + 1)
    us = zeros(T, n, l + 1)

    # Views of the first columns for ease of reference
    u = view(us, :, 1)
    residual = view(rs, :, 1)
    
    # Compute the initial residual rs[:, 1] = b - A * x
    # Avoid computing A * 0.
    if initial_zero
        copy!(residual, b)
    else
        A_mul_B!(residual, A, x)
        @blas! residual -= one(T) * b
        @blas! residual *= -one(T)
        mv_products += 1
    end

    # Apply the left preconditioner
    A_ldiv_B!(Pl, residual)

    # Random shadow residual
    r_shadow = rand(T, n)

    γ = zeros(T, l)
    γ[l] = σ = one(T)

    nrm = norm(residual)
    push!(residuals, nrm)
    iter = 1

    # For the least-squares problem
    M = zeros(T, l + 1, l + 1)
    L = 2 : l + 1

    # Stopping condition based on relative tolerance.
    reltol = nrm * tol

    while nrm > reltol && mv_products < max_mv_products
        σ = -γ[l] * σ
        
        # BiCG part
        for j = 1 : l
            ρ = dot(r_shadow, view(rs, :, j))
            β = ρ / σ
            
            # us[:, 1 : j] .= rs[:, 1 : j] - β * us[:, 1 : j]
            for i = 1 : j
                @blas! view(us, :, i) *= -β
                @blas! view(us, :, i) += one(T) * view(rs, :, i)
            end

            # us[:, j + 1] = Pl \ (A * us[:, j])            
            A_mul_B!(view(us, :, j + 1), A , view(us, :, j))
            A_ldiv_B!(Pl, view(us, :, j + 1))
            mv_products += 1

            σ = dot(r_shadow, view(us, :, j + 1))
            α = ρ / σ

            # rs[:, 1 : j] .= rs[:, 1 : j] - α * us[:, 2 : j + 1]
            for i = 1 : j
                @blas! view(rs, :, i) -= α * view(us, :, i + 1)
            end
            
            # rs[:, j + 1] = Pl \ (A * rs[:, j])
            A_mul_B!(view(rs, :, j + 1), A , view(rs, :, j))
            A_ldiv_B!(Pl, view(rs, :, j + 1))
            mv_products += 1
            
            # x = x + α * us[:, 1]
            @blas! x += α * u
        end

        # MR part: M = rs' * rs
        Ac_mul_B!(M, rs, rs)
        
        # γ = M[L, L] \ M[L, 1] 
        F = lufact!(view(M, L, L))
        A_ldiv_B!(γ, F, view(M, L, 1))

        # This could even be BLAS 3 when combined.
        BLAS.gemv!('N', -one(T), view(us, :, L), γ, one(T), u)
        BLAS.gemv!('N', one(T), view(rs, :, 1 : l), γ, one(T), x)
        BLAS.gemv!('N', -one(T), view(rs, :, L), γ, one(T), residual)
        
        nrm = norm(residual)
        
        visitor(nrm, mv_products)
    end

    x, nrm
end

#################
# Documentation #
#################

let
doc_call = "bicgstab(A, b, l)"
doc!_call = "bicgstab!(x, A, b, l)"

doc_msg = "Solve A*x = b with the BiCGStab(l)"
doc!_msg = "Overwrite `x`.\n\n" * doc_msg

doc_arg = ""
doc!_arg = """`x`: initial guess, overwrite final estimation."""

doc_version = (doc_call, doc_msg, doc_arg)
doc!_version = (doc!_call, doc!_msg, doc!_arg)

docstring = String[]

#Build docs
for (call, msg, arg) in (doc_version, doc!_version) #Start
    push!(docstring, 
"""
$call

$msg

# Arguments

$arg

`A`: linear operator.

`b`: right hand side (vector).

`l::Int = 2`: Number of GMRES steps.

## Keywords

`Pl = Identity()`: left preconditioner of the method.

`tol::Real = sqrt(eps(real(eltype(b))))`: tolerance for stopping condition 
`|r_k| / |r_0| ≤ tol`. Note that:

1. The actual residual is never computed during the iterations; only an 
approximate residual is used.
2. If a preconditioner is given, the stopping condition is based on the 
*preconditioned residual*.

`max_mv_products::Int = min(30, size(A, 1))`: maximum number of matrix
vector products. For BiCGStab this is a less dubious criterion than maximum
number of iterations.

`visitor`: a callback that get's called after each outer iteration, with arguments
`visitor(nrm, mv_products)`. For example use

```julia
logger(nrm, mv_products) = println(mv_products, " ", nrm)
x, residual = bicgstab(A, b, visitor = logger)
```

# Output

`x`: approximated solution.
`residual`: last approximate residual norm

# References
[1] Sleijpen, Gerard LG, and Diederik R. Fokkema. "BiCGstab(l) for 
linear equations involving unsymmetric matrices with complex spectrum." 
Electronic Transactions on Numerical Analysis 1.11 (1993): 2000.
"""
    )
end

@doc docstring[1] -> bicgstabl
@doc docstring[2] -> bicgstabl!
end