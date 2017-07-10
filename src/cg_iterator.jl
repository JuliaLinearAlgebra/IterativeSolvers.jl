import Base: start, next, done

mutable struct CGIterator{matT <: AbstractMatrix, vecT <: AbstractVector, numT <: Real}
    A::matT
    x::vecT
    b::vecT
    r::vecT
    c::vecT
    u::vecT
    reltol::numT
    residual²::numT
    prev_residual²::numT
    maxiter::Int
    mv_products::Int
end

start(it::CGIterator) = 0

done(it::CGIterator, iteration::Int) = iteration ≥ it.maxiter || it.residual² ≤ it.reltol * it.reltol

function next(it::CGIterator, iteration::Int)
    β = it.residual² / it.prev_residual²

    # u := r + βu (almost an axpy)
    @blas! it.u *= β
    @blas! it.u += one(eltype(it.b)) * it.r

    # c = A * u
    A_mul_B!(it.c, it.A, it.u)
    α = it.residual² / dot(it.u, it.c)

    # Improve solution and residual
    @blas! it.x += α * it.u
    @blas! it.r -= α * it.c

    it.prev_residual² = it.residual²
    it.residual² = norm(it.r)^2 # a bit awkward?

    # Return the residual at item and iteration number as state
    √it.residual², iteration + 1
end

@inline cg_iterator(A, b; kwargs...) = cg_iterator!(zeros(b), A, b; initially_zero = true, kwargs...)

function cg_iterator!(x, A, b;
    tol = sqrt(eps(real(eltype(b)))),
    maxiter = min(20, length(b)),
    initially_zero::Bool = false
)
    u = zeros(x)
    r = copy(b)

    # Compute r with an MV-product or not.
    if initially_zero
        c = similar(x)
    else
        c = A * x
        @blas! r -= one(eltype(x)) * c
        mv_products += 1
    end

    # Bookkeeping & stopping criterion
    mv_products = 0
    residual² = norm(r)^2
    prev_residual² = one(residual²)
    reltol = √residual² * tol

    # Return the iterable
    CGIterator(A, x, b,
        r, c, u,
        reltol, residual², prev_residual²,
        maxiter, mv_products
    )
end

function my_cg!(x, A, b; kwargs...)
    it = cg_iterator!(x, A, b; kwargs...)

    for item = it end

    it.x
end