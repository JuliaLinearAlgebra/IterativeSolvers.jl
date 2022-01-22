"""
    LookAheadLanczosDecompOptions

Options for [`LookAheadLanczosDecomp`](@ref).

# Fields
- `max_iter`: Maximum number of iterations
- `max_block_size`: Maximum block size allowed to be constructed
- `log`: Flag determining logging
- `verbose`: Flag determining verbosity
"""
struct LookAheadLanczosDecompOptions
    max_iter::Int
    max_block_size::Int
    log::Bool
    verbose::Bool
end

"""
    LookAheadLanczosDecompLog

Log for [`LookAheadLanczosDecomp`](@ref). In particular, logs the sizes of blocks constructed in the P-Q and V-W sequences.
"""
struct LookAheadLanczosDecompLog
    PQ_block_count::Dict{Int, Int}
    VW_block_count::Dict{Int, Int}
end
LookAheadLanczosDecompLog() = LookAheadLanczosDecompLog(Dict{Int, Int}(), Dict{Int, Int}())

mutable struct LookAheadLanczosDecomp{OpT, OptT, VecT, MatT, ElT, ElRT}
    # Operator
    A::OpT
    At::OptT

    # P-Q sequence
    p::VecT
    q::VecT
    p̂::VecT
    q̂::VecT
    P::MatT
    Q::MatT

    # V-W sequence
    v::VecT
    w::VecT
    ṽ::VecT
    w̃::VecT
    V::MatT
    W::MatT

    # matrix-vector products
    Ap::VecT
    Atq::VecT

    # dot products - note we take tranpose(w)*v, not adjoint(w)*v
    qtAp::ElT
    w̃tṽ::ElT
    wtv::ElT

    # norms
    normp::ElRT
    normq::ElRT
    ρ::ElRT
    ξ::ElRT

    γ::Vector{ElRT}

    D::Matrix{ElT}
    E::Matrix{ElT}
    F::Matrix{ElT}
    F̃lastcol::Vector{ElT}
    G::Matrix{ElT}
    H::Vector{ElT}

    U::UpperTriangular{ElT, Matrix{ElT}}
    L::Matrix{ElT}

    # Indices tracking location in block and sequence
    n::Int
    k::Int
    l::Int
    kstar::Int
    lstar::Int
    mk::Vector{Int}
    nl::Vector{Int}

    # Flag determining if we are in inner block, see [^Freund1994] Alg. 5.2
    innerp::Bool
    innerv::Bool

    # Estimate of norm(A), see [^Freund1993]
    nA::ElRT
    nA_recompute::ElRT

    # Logs and options
    log::LookAheadLanczosDecompLog
    opts::LookAheadLanczosDecompOptions
end

"""
    LookAheadLanczosDecomp(A, v, w; kwargs...)

Provides an iterable which constructs basis vectors for a Krylov subspace generated by `A` given by two initial vectors `v` and `w`. This implementation follows [^Freund1994], where a coupled two-term recurrence is used to generate both a V-W and a P-Q sequence. Following the reference, the Lanczos sequence is generated by `A` and `transpose(A)`.

# Arguments
- `A`: Operator used to construct Krylov subspace
- `v`: Initial vector for Krylov subspace generated by `A`
- `w`: Initial vector for Krylov subspace generated by `transpose(A)`

# Keywords
- `max_iter=size(A, 2)`: Maximum number of iterations 
- `max_block_size=2`: Maximum look-ahead block size to construct. Following [^Freund1994], it is rare for blocks to go beyond size 3. This pre-allocates the block storage for the computation to size `(max_block_size, length(v))`. If a block would be built that exceeds this size, the estimate of `norm(A)` is adjusted to allow the block to close.
- `log=false`: Flag determining whether to log history in a [`LookAheadLanczosDecompLog`](@ref)
- `verbose=false`: Flag determining verbosity of output during iteration

# References
[^Freund1993]:
Freund, R. W., Gutknecht, M. H., & Nachtigal, N. M. (1993). An Implementation of the Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices. SIAM Journal on Scientific Computing, 14(1), 137–158. https://doi.org/10.1137/0914009
[^Freund1994]:
Freund, R. W., & Nachtigal, N. M. (1994). An Implementation of the QMR Method Based on Coupled Two-Term Recurrences. SIAM Journal on Scientific Computing, 15(2), 313–337. https://doi.org/10.1137/0915022
"""
function LookAheadLanczosDecomp(
    A, v, w;
    max_iter=size(A, 2),
    max_block_size=8,
    log=false,
    verbose=false
)
    elT = eltype(v)

    p = similar(v)
    q = similar(v)
    p̂ = similar(v)
    q̂ = similar(v)
    P = similar(v, size(v, 1), 0)
    Q = similar(v, size(v, 1), 0)
    Ap = similar(v)
    Atq = similar(v)
    qtAp = zero(elT)
    normp = zero(real(elT))
    normq = zero(real(elT))

    ṽ = similar(v)
    w̃ = similar(v)
    V = reshape(copy(v), size(v, 1), 1)
    W = reshape(copy(w), size(v, 1), 1)
    w̃tṽ = zero(elT)
    wtv = transpose(w) * v
    ρ = zero(normp)
    ξ = zero(normp)

    γ = Vector{real(elT)}(undef, 1)
    γ[1] = 1.0

    D = Matrix{elT}(undef, 0, 0)
    E = Matrix{elT}(undef, 0, 0)
    G = Matrix{elT}(undef, 0, 0)
    H = Vector{elT}()

    F = Matrix{elT}(undef, 0, 0)
    F̃lastcol = Vector{elT}()

    U = UpperTriangular(Matrix{elT}(undef, 0, 0))
    L = Matrix{elT}(undef, 0, 0)

    n     = 1
    k     = 1
    l     = 1
    kstar = 1
    lstar = 1
    mk    = [1]
    nl    = [1]

    # Equation following 4.5 of [^Freund1993]
    # initial estimate of norm(A)
    nA = max(
        norm(mul!(Ap, A, v)),
        norm(mul!(Atq, transpose(A), w))
    )
    if log
        ld_log = LookAheadLanczosDecompLog(
            Dict(i => 0 for i=1:max_block_size),
            Dict(i => 0 for i=1:max_block_size)
        )
    else
        ld_log = LookAheadLanczosDecompLog()
    end

    return LookAheadLanczosDecomp(
        A, transpose(A),
        p, q, p̂, q̂, P, Q,
        v, w, ṽ, w̃, V, W,
        Ap, Atq,
        qtAp, w̃tṽ, wtv,
        normp, normq, ρ, ξ,
        γ,
        D, E, F, F̃lastcol, G, H,
        U, L,
        n, k, l, kstar, lstar, mk, nl,
        false, false, nA, nA,
        ld_log,
        LookAheadLanczosDecompOptions(
            max_iter,
            max_block_size,
            log,
            verbose
        )
    )
end

isverbose(ld::LookAheadLanczosDecomp) = ld.opts.verbose
islogging(ld::LookAheadLanczosDecomp) = ld.opts.log