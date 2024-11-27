using Printf
import Base: iterate
import LinearAlgebra: UpperTriangular, UpperHessenberg
import BlockDiagonals: BlockDiagonal, blocks
import BlockDiagonals

"""
    LookAheadLanczosDecompOptions

Options for [`IterativeSolvers.LookAheadLanczosDecomp`](@ref).

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

Log for [`IterativeSolvers.LookAheadLanczosDecomp`](@ref). In particular, logs the sizes of blocks constructed in the P-Q and V-W sequences.
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
    P::LimitedMemoryMatrix{ElT, MatT}
    Q::LimitedMemoryMatrix{ElT, MatT}

    # V-W sequence
    v::VecT
    w::VecT
    ṽ::VecT
    w̃::VecT
    V::LimitedMemoryMatrix{ElT, MatT}
    W::LimitedMemoryMatrix{ElT, MatT}

    # matrix-vector products
    Ap::VecT
    Atq::VecT

    # dot products - note we take tranpose(w)*v, not adjoint(w)*v
    qtAp::ElT
    w̃tṽ::ElT
    wtv::ElT

    # norms
    # We store normp and normq to aid in checks Eq. 4.6 and 4.7
    normp::Vector{ElRT}
    normq::Vector{ElRT}
    ρ::ElRT
    ξ::ElRT

    γ::Vector{ElRT}

    # Eq. 2.13
    D::BlockDiagonal{ElT, Matrix{ElT}}
    # Eq. 3.14
    E::BlockDiagonal{ElT, Matrix{ElT}}
    # Defined after Eq. 5.1
    Flastcol::Vector{ElT} # size n
    Flastrow::Vector{ElT} # size n-1
    F̃lastcol::Vector{ElT}
    # Eq. 5.1
    G::Vector{ElT}
    # Eq. 3.11
    H::Vector{ElT}

    # Eq. 3.9
    # need to keep previous columns of U for G checks
    U::LimitedMemoryUpperTriangular{ElT, Matrix{ElT}}
    # need to keep previous columns of L for H checks
    L::LimitedMemoryUpperHessenberg{ElT, Matrix{ElT}}
    
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

    # Logs and options
    log::LookAheadLanczosDecompLog
    opts::LookAheadLanczosDecompOptions
end

"""
    LookAheadLanczosDecomp(A, v, w; kwargs...)

Provides an iterable which constructs basis vectors for a Krylov subspace generated by `A` given by two initial vectors `v` and `w` [^Freund1993]. This implementation follows [^Freund1994], where a coupled two-term recurrence is used to generate both a V-W and a P-Q sequence. Following the reference, the Lanczos sequence is generated by `A` and `transpose(A)`.

# Arguments
- `A`: Operator used to construct Krylov subspace
- `v`: Initial vector for Krylov subspace generated by `A`
- `w`: Initial vector for Krylov subspace generated by `transpose(A)`

# Keywords
- `vw_normalized=false`: Flag if `v` and `w` passed to function are already normalized. If `false`, normalized `v` and `w` in-place.
- `max_iter=size(A, 2)`: Maximum number of iterations 
- `max_block_size=4`: Maximum look-ahead block size to construct. Following [^Freund1994], it is rare for blocks to go beyond size 3. This pre-allocates the block storage for the computation to size `(max_block_size, length(v))`. If a block would be built that exceeds this size, the estimate of `norm(A)` is adjusted to allow the block to close.
- `max_memory=4`: Maximum memory to store the sequence vectors. This may be greater than the block size.
- `log=false`: Flag determining whether to log history in a [`IterativeSolvers.LookAheadLanczosDecompLog`](@ref)
- `verbose=false`: Flag determining verbosity of output during iteration

# References
[^Freund1993]:
Freund, R. W., Gutknecht, M. H., & Nachtigal, N. M. (1993). An Implementation of the Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices. SIAM Journal on Scientific Computing, 14(1), 137–158. https://doi.org/10.1137/0914009

[^Freund1994]:
Freund, R. W., & Nachtigal, N. M. (1994). An Implementation of the QMR Method Based on Coupled Two-Term Recurrences. SIAM Journal on Scientific Computing, 15(2), 313–337. https://doi.org/10.1137/0915022
"""
function LookAheadLanczosDecomp(
    A, v, w;
    vw_normalized=false,
    max_iter=size(A, 2),
    max_block_size=4,
    max_memory=4,
    log=false,
    verbose=false
)
    elT = eltype(v)

    # Alg 5.2.0
    if !vw_normalized
        normalize!(v)
        normalize!(w)
    end

    p = similar(v)
    q = similar(v)
    p̂ = similar(v)
    q̂ = similar(v)
    P = LimitedMemoryMatrix(similar(v, size(v, 1), 0), max_memory)
    Q = LimitedMemoryMatrix(similar(v, size(v, 1), 0), max_memory)
    Ap = similar(v)
    Atq = similar(v)
    qtAp = zero(elT)
    normp = Vector{real(elT)}()
    normq = Vector{real(elT)}()

    ṽ = similar(v)
    w̃ = similar(v)
    V = LimitedMemoryMatrix(copy(v), max_memory)
    W = LimitedMemoryMatrix(copy(w), max_memory)
    w̃tṽ = zero(elT)
    
    wtv = transpose(w) * v

    ρ = zero(real(elT))
    ξ = zero(real(elT))

    γ = Vector{real(elT)}(undef, 1)
    γ[1] = 1.0

    D = BlockDiagonal{elT, Matrix{elT}}(Vector{Matrix{elT}}())
    E = BlockDiagonal{elT, Matrix{elT}}(Vector{Matrix{elT}}())
    G = Vector{elT}()
    H = Vector{elT}()

    Flastcol = Vector{elT}()
    Flastrow = Vector{elT}()
    F̃lastcol = Vector{elT}()

    U = LimitedMemoryUpperTriangular{elT, Matrix{elT}}(max_memory)
    L = LimitedMemoryUpperHessenberg{elT, Matrix{elT}}(max_memory)

    # Alg 5.2.0
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
        D, E, Flastcol, Flastrow, F̃lastcol, G, H,
        U, L,
        n, k, l, kstar, lstar, mk, nl,
        false, false, nA,
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

# Note: Equation references are w.r.t [^Freund1994] unless otherwise stated
# Definitions and Eq. 3.27, 3.28
_PQ_block_size(ld) = max(1, ld.n - ld.mk[ld.k])
# Definitions and Eq. 2.20, 2.21
_VW_block_size(ld) = ld.n+1 - ld.nl[ld.l]
_VW_prev_block_size(ld) = ld.nl[ld.l] - ld.nl[max(1, ld.l-1)]
_is_block_small(ld, n) = n < ld.opts.max_block_size

"""
    _grow_last_block!(A, Bcol, Brow, Bcorner)

Grows the last block in-place in `A` by appending the column `Bcol`, the row `Brow`, and the corner element `Bcorner`. `Bcol` and `Brow` are automatically truncated to match the size of the grown block
"""
function _grow_last_block!(A::BlockDiagonal{T, TM}, Bcol, Brow, Bcorner) where {T, TM}
    n = BlockDiagonals.nblocks(A)
    b = BlockDiagonals.blocks(A)
    s = size(last(b), 1)
    b[n] = TM([
        b[n] Bcol[end-s+1:end]
        Brow[1:1, end-s+1:end] Bcorner
    ])
    return A
end

"""
    _start_new_block!(A, B)

Appends a new block to the end of `A` with `B`
"""
function _start_new_block!(A::BlockDiagonal{T, TM}, B) where {T, TM}
    push!(BlockDiagonals.blocks(A), TM(fill(only(B), 1, 1)))
    return A
end

start(::LookAheadLanczosDecomp) = 1
done(ld::LookAheadLanczosDecomp, iteration::Int) = iteration ≥ ld.opts.max_iter
function iterate(ld::LookAheadLanczosDecomp, n::Int=start(ld))
    isverbose(ld) && n==1 && @printf("%7s\t%7s\t%7s\t%7s\t%7s\t%7s\n", "n", "γ", "|p|", "|q|", "ρ", "ξ")
    # Alg. 5.2.1 - 5.2.11
    _update_PQ_sequence!(ld)
    # Alg. 5.2.12
    if ld.normp[end] < eps(eltype(ld.normp)) || ld.normq[end] < eps(eltype(ld.normq))
        return nothing
    end
    # Alg. 5.2.13 - 5.2.25
    ld, finished = _update_VW_sequence!(ld)
    isverbose(ld) && @printf("%3d\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\n", n, ld.γ[end], ld.normp[end], ld.normq[end], ld.ρ[end], ld.ξ[end])
    if finished return nothing end

    # if on final iteration, finish
    if done(ld, n) return nothing end

    # Loop condition
    ld.n += 1
    return ld, ld.n
end

function _update_PQ_sequence!(ld)
    # check if we can make an inner block
    inner_ok = _is_block_small(ld, _PQ_block_size(ld))
    # Alg. 5.2.1
    _update_D!(ld)
    # Alg. 5.2.2
    _update_kstar!(ld)
    # Alg. 5.2.3
    _update_Flastrow!(ld)
    # Alg. 5.2.4
    ld.innerp = inner_ok && !isempty(blocks(ld.E)) && _is_singular(last(blocks(ld.E)))
    # Alg. 5.2.5
    _update_U!(ld, ld.innerp)
    # Alg. 5.2.6
    _update_p̂q̂_common!(ld)
    # Condition from Alg. 5.2.6
    if !ld.innerp
        # Alg. 5.2.7
        _update_Gnm1!(ld)
        ld.innerp = inner_ok && _check_G(ld)
        # Condition from Alg. 5.2.7 - we are confident we can perform a regular step
        if !ld.innerp
            # Alg. 5.2.8
            _update_pq_regular!(ld)
            _matvec_pq!(ld)
            # Alg. 5.2.9
            _update_Gn!(ld)
            ld.innerp = inner_ok && _check_G(ld)
            # Condition from Alg. 5.2.9 - One last chance to bail and create an inner step
            if !ld.innerp
                if islogging(ld)
                    ld.log.PQ_block_count[_PQ_block_size(ld)] += 1
                end
                if !isone(ld.n)
                    # Alg. 5.2.10
                    push!(ld.mk, ld.n)
                    ld.k += 1
                end
            else
                # Alg. 5.2.11
                # re-calculate matrix-vector product and append new result to P, Q seq
                isverbose(ld) && @info "Inner P-Q construction, second G check"
                _update_pq_inner!(ld)
                _matvec_pq!(ld)
            end
        else
            # Alg. 5.2.11
            isverbose(ld) && @info "Inner P-Q construction, first G check"
            _update_pq_inner!(ld)
            _matvec_pq!(ld)
        end
    else
        # Alg. 5.2.11
        isverbose(ld) && @info "Inner P-Q construction, singular E check"
        _update_pq_inner!(ld)
        _matvec_pq!(ld)
    end
    _append_PQ!(ld)
    return ld
end

function _update_VW_sequence!(ld)
    # check if we can make an inner block
    inner_ok = _is_block_small(ld, _VW_block_size(ld))
    # Alg. 5.2.13
    mul!(ld.Atq, ld.At, ld.q)
    # Alg. 5.2.14
    _update_E!(ld)
    # Alg. 5.2.15
    _update_lstar!(ld)
    # Alg. 5.2.16
    _update_Flastcol!(ld)
    # Alg. 5.2.17
    ld.innerv = inner_ok && _is_singular(last(blocks(ld.D)))
    # Alg. 5.2.18
    _update_L!(ld, ld.innerv)
    # Alg. 5.2.19
    _update_v̂ŵ_common!(ld)
    # Condition from # Alg. 5.2.19
    if !ld.innerv
        # Alg. 5.2.20
        _update_Hn!(ld)
        ld.innerv = inner_ok && _check_H(ld)
        if !ld.innerv
            # Alg. 5.2.21
            _update_vw_regular!(ld)
            # also updates γ, ω̃tṽ
            ld, terminate_early = _update_ρξ!(ld)
            if terminate_early return ld, true end
            # Alg. 5.2.22
            _update_Hnp1!(ld)
            ld.innerv = inner_ok && _check_H(ld)
            # Condition from Alg. 5.2.22
            if !ld.innerv
                if islogging(ld)
                    ld.log.VW_block_count[_VW_block_size(ld)] += 1
                end
                # Alg. 5.2.23
                push!(ld.nl, ld.n+1)
                ld.l += 1
            else
                # Alg. 5.2.24
                isverbose(ld) && @info "Inner V-W construction, second H check"
                _update_vw_inner!(ld)
                ld.innerv = true
                ld, terminate_early = _update_ρξ!(ld, true)
                if terminate_early return ld, true end
            end
        else
            # Alg. 5.2.24
            isverbose(ld) && @info "Inner V-W construction, first H check"
            _update_vw_inner!(ld)
            ld, terminate_early = _update_ρξ!(ld)
            if terminate_early return ld, true end
        end
    else
        # Alg. 5.2.24
        isverbose(ld) && @info "Inner V-W construction, singular D check"
        _update_vw_inner!(ld)
        ld, terminate_early = _update_ρξ!(ld)
        if terminate_early return ld, true end
    end
    # Alg. 5.2.25
    _update_vw!(ld)
    return ld, false
end

function _update_D!(ld)
    # Alg. 5.2.1
    # Eq. 5.2:
    # F[n-1] = Wt[n-1]V[n]L[n-1] = D[n-1]L[1:n-1, 1:n-1] + l[n, n-1]D[1:n-1, n][0 ... 0 1]
    # => D[1:end-1, end] = (F[:, end] - (D_prev L[1:end-1, 1:end]))[:, end] / ρ
    # Eq. 3.15, (D Γ)ᵀ = (D Γ)
    # D[n, n] = wtv

    if isone(ld.n) || _VW_block_size(ld) == 1
        _start_new_block!(ld.D, ld.wtv)
    else
        Dlastcol = (ld.Flastcol - (ld.D * ld.L[1:end-1, end])) / ld.ρ
        Dlastrow = transpose(Dlastcol * ld.γ[end] ./ ld.γ[1:end-1])
        _grow_last_block!(ld.D, Dlastcol, Dlastrow, ld.wtv)
    end
    return ld
end

function _update_kstar!(ld)
    # Alg. 5.2.2
    # Eq. 3.31
    # Returns the largest block index in the P-Q sequence which starts at the same index
    # as the current V-W block with index l
    k, mk, l, nl = ld.k, ld.mk, ld.l, ld.nl
    ld.kstar = max(1, min(k, searchsortedlast(mk, nl[l]-1)))
    return ld
end

function _update_Flastrow!(ld)
    # Alg. 5.2.3
    # Auxiliary matrix F = WtAP, F̃ = QtAV
    # Lemma 5.1: FΓ = (F̃Γ)ᵀ
    # Note: D is at D_n, L is at L_{n-1}
    # Eq. 5.2 (w/ indices advanced): 
    # F_{n} = D_{n}L[1:n, 1:n] + l[n+1, n]D_{n}[1:n, n+1][0 ... 0 1]
    if !isone(ld.n) # We only need to do this if we are constructing a block
        ld.Flastrow = reshape(ld.D[end:end, :] * ld.L, :)
        ld.F̃lastcol = ld.Flastrow .* ld.γ[1:end-1] ./ ld.γ[end]
    end
end

function _update_U!(ld, innerp)
    # Alg. 5.2.5
    # U is upper triangular matrix in decomposition of recurrence relation for P-Q sequence
    # updates last column of U
    n, mk, k, kstar = ld.n, ld.mk, ld.k, ld.kstar
    uvec = fill(0, n)
    uvec[end] = 1
    _grow_hcat!(ld.U, uvec)

    for i = kstar:k-1
        block_start = mk[i]
        block_end = mk[i+1]-1
        ld.U[block_start:block_end, end] .= blocks(ld.E)[i] \ ld.F̃lastcol[block_start:block_end]
    end
    if !innerp && !isone(n)
        ld.U[mk[k]:end-1, end] .= blocks(ld.E)[end] \ ld.F̃lastcol[mk[k]:end]
    end
    return ld
end

function _update_p̂q̂_common!(ld)
    # Alg. 5.2.6
    mk, k, kstar = ld.mk, ld.k, ld.kstar
    copyto!(ld.p̂, ld.v)
    copyto!(ld.q̂, ld.w)
    for i = mk[kstar]:mk[k]-1 # TODO: OPTIMIZE gemv! (or 5-arg mul!)
        if ld.U[i, end] != 0
            axpy!(-ld.U[i, end], view(ld.P, :, i), ld.p̂)
            axpy!(-ld.U[i, end] * ld.γ[end] / ld.γ[i], view(ld.Q, :, i), ld.q̂)
        end
    end
end
function _update_Gnm1!(ld)
    # Alg. 5.2.7
    # G_{n-1} = U_n L_{n-1}
    # U is currently U_n and L is currently L_{n-1}
    n, mk, k = ld.n, ld.mk, ld.k
    ld.G = fill(0.0, n-1)
    if !isone(n)
        ld.G[ld.mk[ld.k]:end] .= ld.U[mk[k]:end-1, mk[k]:end] * ld.L[mk[k]:end, end]
    end
    return ld
end

function _update_Gn!(ld)
    # G_{n-1} = U_n L_{n-1}
    # G[m_k:n-1, n] = E_{k} \ Q_{k}ᵀ A Ap_n
    # Q_{n-1}ᵀ A Ap_n = γ_{n-1}/γ_n ρ_{n} [0 ⋯ 0 q_nᵀAp_n]ᵀ
    n, mk, k = ld.n, ld.mk, ld.k
    if !isone(ld.n)
        qtAp = fill(zero(eltype(ld.G)), length(ld.G))
        qtAp[end] = ld.qtAp
        ld.G = (ld.E \ qtAp * ld.ρ * ld.γ[n-1] / ld.γ[n])
    end

    return ld
end

function _check_G(ld)
    # Eq. 4.6, 4.7
    n = ld.n
    if n <= 2 return false end
    return !(
        ld.nA * ld.normp[end] ≥ sum(abs(ld.G[i]) * ld.normp[i] for i in 1:length(ld.normp)-1) &&
        ld.nA * ld.normq[end] ≥ sum(ld.γ[n-1]/ld.γ[i] * abs(ld.G[i]) * ld.normq[i] for i in 1:length(ld.normq)-1)
    )
    return false
end

function _update_pq_regular!(ld)
    # 5.2.8
    n, mk, k, kstar = ld.n, ld.mk, ld.k, ld.kstar
    copyto!(ld.p, ld.p̂)
    copyto!(ld.q, ld.q̂)
    for i = mk[k]:n-1 # TODO: OPTIMIZE gemv! (or 5-arg mul!)
        if ld.U[i, end] != 0
            axpy!(-ld.U[i, end], view(ld.P, :, i), ld.p)
            axpy!(-ld.U[i, end] * ld.γ[n] / ld.γ[i], view(ld.Q, :, i), ld.q)
        end
    end
    return ld
end

function _update_pq_inner!(ld)
    # 5.2.11
    n, mk, k, kstar = ld.n, ld.mk, ld.k, ld.kstar
    copyto!(ld.p, ld.p̂)
    copyto!(ld.q, ld.q̂)
    for i = mk[k]:n-1 # TODO: OPTIMIZE gemv!
        u = _u(i, n, mk[k])
        ld.U[i, end] = u
        if u != 0
            axpy!(-u, view(ld.P, :, i), ld.p)
            axpy!(-u * ld.γ[n] / ld.γ[i], view(ld.Q, :, i), ld.q)
        end
    end
    return ld
end

function _matvec_pq!(ld)
    # Common part of Alg. 5.2.8, Alg. 5.2.11
    mul!(ld.Ap, ld.A, ld.p)
    ld.qtAp = transpose(ld.q) * ld.Ap
    return ld
end

function _append_PQ!(ld)
    # adds p, q vector sequence to P, Q
    hcat!(ld.P, ld.p)
    hcat!(ld.Q, ld.q)
    push!(ld.normp, norm(ld.p))
    push!(ld.normq, norm(ld.q))
end

# 5.2.4, 5.2.17
# Freund, R. W., Gutknecht, M. H., & Nachtigal, N. M. (1993). An Implementation of the Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices Part I. SIAM Journal on Scientific Computing, 14(1), 137–158. https://doi.org/10.1137/0914009
# suggests eps()^(1/3)
_is_singular(A) = !isempty(A) && minimum(svdvals(A)) < eps(real(eltype(A)))^(1/3)

function _update_E!(ld)
    # E = QtAP
    # F = ΓUtinv(Γ)E
    # 5.2.14
    n = ld.n
    if isone(ld.n) || (ld.n == ld.mk[end])
        _start_new_block!(ld.E, ld.qtAp)
    else
        ΓUtinvΓ = ld.γ .* transpose(ld.U) ./ transpose(ld.γ)
        Elastrow = (ΓUtinvΓ[end, end] \ reshape(ld.Flastrow, 1, :) - ΓUtinvΓ[end:end, 1:end-1]*ld.E)
        Elastcol = (transpose(Elastrow) .* ld.γ[1:n-1] ./ ld.γ[n])
        _grow_last_block!(ld.E, Elastcol, Elastrow, ld.qtAp)
    end
    return ld
end

function _update_lstar!(ld)
    # 5.2.15
    # Returns the largest block index in the V-W sequence which starts at the same index
    # as the current P-Q block with index k
    k, mk, l, nl = ld.k, ld.mk, ld.l, ld.nl
    ld.lstar = min(l, searchsortedlast(nl, mk[k]))
end

function _update_Flastcol!(ld)
    # 5.2.16
    n = ld.n
    ΓUtinvΓ = ld.γ .* transpose(ld.U) ./ transpose(ld.γ)
    # length n, ld.F_lastrow of length n-1
    if isone(n)
        ld.Flastcol = fill(ΓUtinvΓ[end, end] * ld.E[end, end], 1)
    else
        ld.Flastcol = ΓUtinvΓ * ld.E[:, end]
    end
    return ld
end

function _update_L!(ld, innerv)
    # Alg. 5.2.18
    n, l, lstar, nl = ld.n, ld.l, ld.lstar, ld.nl
    Llastcol = fill(zero(eltype(ld.L)), n)
    for i = lstar:l-1
        block_start = nl[i]
        block_end = nl[i+1]-1
        Llastcol[block_start:block_end] .= blocks(ld.D)[i] \ ld.Flastcol[block_start:block_end]
    end
    if !innerv
        Llastcol[nl[l]:end] .= blocks(ld.D)[end] \ ld.Flastcol[nl[l]:end]
    end
    lvec = [Llastcol; 0]
    _grow_hcat!(ld.L, lvec)
    return ld
end

function _update_v̂ŵ_common!(ld)
    # 5.2.6
    l, lstar, nl, n = ld.l, ld.lstar, ld.nl, ld.n

    for i = nl[lstar]:nl[l]-1 # TODO: OPTIMIZE gemv! (or 5-arg mul!)
        if ld.L[i, n] != 0
            axpy!(-ld.L[i, n], view(ld.V, :, i), ld.Ap)
            axpy!(-ld.L[i, n] * ld.γ[n] / ld.γ[i], view(ld.W, :, i), ld.Atq)
        end
    end
    return ld
end

function _update_Hn!(ld)
    # Alg. 5.2.20
    n, l, lstar, nl = ld.n, ld.l, ld.lstar, ld.nl
    ld.H = fill(0.0, n)
    if !isone(ld.n)
        mul!(ld.H[nl[l]:end], ld.L[nl[l]:end-1, nl[l]:end], ld.U[nl[l]:end, end])
    end
    return ld
end

function _update_Hnp1!(ld)
    # Alg. 5.2.22
    # note that γ, ρ, w̃tṽ, and ξ are at n+1
    n = ld.n
    wv = fill(zero(eltype(ld.H)), 1, length(ld.H))
    wv[end] = ld.w̃tṽ / (ld.ρ * ld.ξ)
    ld.H .= reshape(ld.D \ transpose(wv) * ld.ρ * ld.γ[n] / ld.γ[n+1], length(wv))
    return ld
end

function _check_H(ld)
    # Eq. 3.11: H = LU
    # Eq. 4.3, 4.4
    n = ld.n
    return !(
        ld.nA ≥ sum(abs(ld.H[i, end]) for i in 1:n) &&
        ld.nA ≥ sum(abs(ld.H[i, end]) * ld.γ[n] / ld.γ[i] for i in 1:n)
    )
end

function _update_vw_regular!(ld)
    # 5.2.8
    n, l, lstar, nl, k, mk = ld.n, ld.l, ld.lstar, ld.nl, ld.k, ld.mk
    copyto!(ld.ṽ, ld.Ap)
    copyto!(ld.w̃, ld.Atq)
    for i = nl[l]:n # TODO: OPTIMIZE gemv! (or 5-arg mul!)
        if ld.L[i, end] != 0
            axpy!(-ld.L[i, end], view(ld.V, :, i), ld.ṽ)
            axpy!(-ld.L[i, end] * ld.γ[n] / ld.γ[i], view(ld.W, :, i), ld.w̃)
        end
    end
    return ld
end

function _update_vw_inner!(ld)
    # 5.2.11
    n, l, k, lstar, nl, mk = ld.n, ld.l, ld.k, ld.lstar, ld.nl, ld.mk
    copyto!(ld.ṽ, ld.Ap)
    copyto!(ld.w̃, ld.Atq)
    for i = nl[l]:n # TODO: OPTIMIZE gemv!
        ll = _l(i, n, nl[l])
        ld.L[i, end] = ll
        if ll != 0
            axpy!(-_l(i, n, nl[l]), view(ld.V, :, i), ld.ṽ)
            axpy!(-_l(i, n, nl[l]) * ld.γ[n] / ld.γ[i], view(ld.W, :, i), ld.w̃)
        end
    end
    return ld
end

function _update_ρξ!(ld, retry=false)
    # Common part of Alg. 5.2.21 and Alg. 5.2.24
    # if retry, then this means we have already added data to the vectors, but our
    # inner block check failed, so we overwrite what he have. This is the case if
    # the check at Alg. 5.2.22 fails
    n = ld.n
    ld.ρ = norm(ld.ṽ)
    ld.L[end, end] = ld.ρ
    ld.ξ = norm(ld.w̃)
    terminate_early = false
    if ld.ρ < eps(typeof(ld.ρ)) || ld.ξ < eps(typeof(ld.ξ))
        terminate_early = true
    else
        if retry
            ld.γ[n+1] = ld.γ[n]*ld.ρ/ld.ξ
        else
            push!(ld.γ, ld.γ[n]*ld.ρ/ld.ξ)
        end
        ld.w̃tṽ = transpose(ld.w̃) * ld.ṽ
    end
    return ld, terminate_early
end

function _update_vw!(ld)
    # 5.2.25
    ld.v = ld.ṽ / ld.ρ
    ld.w = ld.w̃ / ld.ξ
    ld.wtv = ld.w̃tṽ / (ld.ρ * ld.ξ)
    hcat!(ld.V, ld.v)
    hcat!(ld.W, ld.w)
    return ld
end

function _u(i, n, mk)
    # inner recurrence coefficients
    if i == n-1
        return 1.0
    elseif i == n-2 && mk <= n-2
        return 1.0
    else
        return 0.0
    end
end
function _l(i, n, nl)
    # inner recurrence coefficients
    if i == n
        return 1.0
    elseif i == n-1 && nl <= n-1
        return 1.0
    else
        return 0.0
    end
end