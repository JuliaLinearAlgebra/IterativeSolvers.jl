#=
The code below was derived from the scipy implementation of the LOBPCG algorithm in https://github.com/scipy/scipy/blob/v1.1.0/scipy/sparse/linalg/eigen/lobpcg/lobpcg.py#L109-L568. 

Since the link above mentions the license is BSD license, the notice for the BSD license 2.0 is hereby given below giving credit to the authors of the Python implementation.

Copyright (c) 2018, Robert Cimrman, Andrew Knyazev
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * The names of the contributors may not be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#

export lobpcg, lobpcg!, LOBPCGIterator

struct LOBPCGState{TR,TL}
    iteration::Int
    residual_norms::TR
    ritz_values::TL
end
function Base.show(io::IO, t::LOBPCGState)
    @printf io "%8d    %14e\n" t.iteration maximum(t.residual_norms)
    return
end

const LOBPCGTrace{TR,TL} = Vector{LOBPCGState{TR,TL}}
function Base.show(io::IO, tr::LOBPCGTrace)
    @printf io "Iteration    Maximum residual norm \n"
    @printf io "---------    ---------------------\n"
    for state in tr
        show(io, state)
    end
    return
end

struct LOBPCGResults{T, TL, TR, TX, TI, TTrace <: LOBPCGTrace}
    λ::TL
    X::TX
    iterations::Int
    residual_norms::TR
    tolerance::T
    converged::Bool
    maxiter::TI
    trace::TTrace
end
function Base.show(io::IO, r::LOBPCGResults)
    first_two(fr) = [x for (i, x) in enumerate(fr)][1:2]

    @printf io "Results of LOBPCG Algorithm\n"
    @printf io " * Algorithm: LOBPCG - CholQR\n"

    if length(join(r.λ, ",")) < 40
        @printf io " * λ: [%s]\n" join(r.λ, ",")
    else
        @printf io " * λ: [%s, ...]\n" join(first_two(r.λ), ",")
    end

    if length(join(r.residual_norms, ",")) < 40
        @printf io " * Residual norm(s): [%s]\n" join(r.residual_norms, ",")
    else
        @printf io " * Residual norm(s): [%s, ...]\n" join(first_two(r.residual_norms), ",")
    end
    @printf io " * Convergence\n"
    @printf io "   * Iterations: %s\n" r.iterations
    @printf io "   * Converged: %s\n" r.converged
    @printf io "   * Iterations limit: %s\n" r.maxiter
    @printf io "   * Reached Maximum Number of Iterations: %s\n" r.iterations >= r.maxiter

    return
end

struct Blocks{Generalized, T, TA<:AbstractArray{T}}
    block::TA # X, R or P
    A_block::TA # AX, AR or AP
    B_block::TA # BX, BR or BP
end
Blocks(X, AX) = Blocks{false, eltype(X), typeof(X)}(X, AX, X)
Blocks(X, AX, BX) = Blocks{true, eltype(X), typeof(X)}(X, AX, BX)
function A_mul_X!(b::Blocks, A)
    A_mul_B!(b.A_block, A, b.block)
    return
end
function A_mul_X!(b::Blocks, A, n)
    A_mul_B!(view(b.A_block, :, 1:n), A, view(b.block, :, 1:n))
    return
end
function B_mul_X!(b::Blocks{true}, B)
    A_mul_B!(b.B_block, B, b.block)
    return
end
function B_mul_X!(b::Blocks{true}, B, n)
    A_mul_B!(view(b.B_block, :, 1:n), B, view(b.block, :, 1:n))
    return
end
function B_mul_X!(b::Blocks{false}, B, n = 0)
    return
end

struct Constraint{T, TA<:AbstractArray{T}, TC}
    Y::TA
    BY::TA
    gram_chol::TC
    gramYBV::TA # to be used in view
    tmp::TA # to be used in view
end
function Constraint(::Void, B, X)
    return Constraint{Void, Matrix{Void}, Void}(Matrix{Void}(0,0), Matrix{Void}(0,0), nothing, Matrix{Void}(0,0), Matrix{Void}(0,0))
end
function Constraint(Y, B, X)
    T = eltype(X)
    if B isa Void
        B = Y
    else
        BY = similar(B)
        A_mul_B!(BY, B, Y)
    end
    gramYBY = Ac_mul_B(Y, BY)
    realdiag!(gramYBY)
    gramYBY_chol = cholfact!(Hermitian(gramYBY))
    gramYBV = zeros(T, size(Y, 2), size(X, 2))
    tmp = similar(gramYBV)

    return Constraint(Y, BY, gramYBY_chol, gramYBV, tmp)
end

function (constr!::Constraint{Void})(X)
    nothing
end

function (constr!::Constraint)(X)
    sizeX = size(X, 2)
    sizeY = size(constr!.Y, 2)
    gramYBV_view = view(constr!.gramYBV, 1:sizeY, 1:sizeX)
    Ac_mul_B!(gramYBV_view, constr!.BY, X)
    tmp_view = view(constr!.tmp, 1:sizeY, 1:sizeX)
    A_ldiv_B!(tmp_view, gram_chol, gramYBV_view)
    A_mul_B!(X, constr!.Y, tmp_view)

    nothing
end

struct RPreconditioner{TM, T, TA<:AbstractArray{T}}
    M::TM
    buffer::TA
    RPreconditioner{TM, T, TA}(M, X) where {TM, T, TA<:AbstractArray{T}} = new(M, similar(X))
end
RPreconditioner(M, X) = RPreconditioner{typeof(M), eltype(X), typeof(X)}(M, X)

function (precond!::RPreconditioner{Void})(X)
    nothing
end
function (precond!::RPreconditioner)(X)
    bs = size(X, 2)
    A_ldiv_B!(view(precond!.buffer, :, 1:bs), precond!.M, X)
    # Just returning buffer would be cheaper but struct at call site must be mutable
    X .= view(precond!.buffer, :, 1:bs)
    nothing
end

struct BlockGram{Generalized, TA}
    XAX::TA
    XAR::TA
    XAP::TA
    RAR::TA
    RAP::TA
    PAP::TA
end
function BlockGram(XBlocks::Blocks{Generalized, T}) where {Generalized, T}
    sizeX = size(XBlocks.block, 2)
    XAX = zeros(T, sizeX, sizeX)
    XAP = zeros(T, sizeX, sizeX)
    XAR = zeros(T, sizeX, sizeX)
    RAR = zeros(T, sizeX, sizeX)
    RAP = zeros(T, sizeX, sizeX)
    PAP = zeros(T, sizeX, sizeX)
    return BlockGram{Generalized, Matrix{T}}(XAX, XAR, XAP, RAR, RAP, PAP)
end
XAX!(BlockGram, XBlocks) = Ac_mul_B!(BlockGram.XAX, XBlocks.block, XBlocks.A_block)
XAP!(BlockGram, XBlocks, PBlocks, n) = Ac_mul_B!(view(BlockGram.XAP, :, 1:n), XBlocks.block, view(PBlocks.A_block, :, 1:n))
XAR!(BlockGram, XBlocks, RBlocks, n) = Ac_mul_B!(view(BlockGram.XAR, :, 1:n), XBlocks.block, view(RBlocks.A_block, :, 1:n))
RAR!(BlockGram, RBlocks, n) = Ac_mul_B!(view(BlockGram.RAR, 1:n, 1:n), view(RBlocks.block, :, 1:n), view(RBlocks.A_block, :, 1:n))
RAP!(BlockGram, RBlocks, PBlocks, n) = Ac_mul_B!(view(BlockGram.RAP, 1:n, 1:n), view(RBlocks.A_block, :, 1:n), view(PBlocks.block, :, 1:n))
PAP!(BlockGram, PBlocks, n) = Ac_mul_B!(view(BlockGram.PAP, 1:n, 1:n), view(PBlocks.block, :, 1:n), view(PBlocks.A_block, :, 1:n))
XBP!(BlockGram, XBlocks, PBlocks, n) = Ac_mul_B!(view(BlockGram.XAP, :, 1:n), XBlocks.block, view(PBlocks.B_block, :, 1:n))
XBR!(BlockGram, XBlocks, RBlocks, n) = Ac_mul_B!(view(BlockGram.XAR, :, 1:n), XBlocks.block, view(RBlocks.B_block, :, 1:n))
RBP!(BlockGram, RBlocks, PBlocks, n) = Ac_mul_B!(view(BlockGram.RAP, 1:n, 1:n), view(RBlocks.B_block, :, 1:n), view(PBlocks.block, :, 1:n))
XBX!(BlockGram, XBlocks) = Ac_mul_B!(BlockGram.XAX, XBlocks.block, XBlocks.B_block)
RBR!(BlockGram, RBlocks, n) = Ac_mul_B!(view(BlockGram.RAR, 1:n, 1:n), view(RBlocks.block, :, 1:n), view(RBlocks.B_block, :, 1:n))
PBP!(BlockGram, PBlocks, n) = Ac_mul_B!(view(BlockGram.PAP, 1:n, 1:n), view(PBlocks.block, :, 1:n), view(PBlocks.B_block, :, 1:n))

function I!(G, xr)
    for j in xr, i in xr
        G[i, j] = ifelse(i==j, 1, 0)
    end
    return
end

function (g::BlockGram)(gram, lambda, n1::Int, n2::Int, n3::Int)
    xr = 1:n1
    rr = n1+1:n1+n2
    pr = n1+n2+1:n1+n2+n3
    if n1 > 0
        #gram[xr, xr] .= view(g.XAX, 1:n1, 1:n1)
        gram[xr, xr] .= Diagonal(view(lambda, 1:n1))
    end
    if n2 > 0
        gram[rr, rr] .= view(g.RAR, 1:n2, 1:n2)
        gram[xr, rr] .= view(g.XAR, 1:n1, 1:n2)
        transpose!(view(gram, rr, xr), view(g.XAR, 1:n1, 1:n2))
    end
    if n3 > 0
        gram[pr, pr] .= view(g.PAP, 1:n3, 1:n3)
        gram[rr, pr] .= view(g.RAP, 1:n2, 1:n3)
        gram[xr, pr] .= view(g.XAP, 1:n1, 1:n3)
        transpose!(view(gram, pr, rr), view(g.RAP, 1:n2, 1:n3))
        transpose!(view(gram, pr, xr), view(g.XAP, 1:n1, 1:n3))
    end
    return 
end
function (g::BlockGram)(gram, n1::Int, n2::Int, n3::Int, normalized::Bool=true)
    xr = 1:n1
    rr = n1+1:n1+n2
    pr = n1+n2+1:n1+n2+n3
    if n1 > 0
        if normalized
            I!(gram, xr)
        else
            gram[xr, xr] .= view(g.XAX, 1:n1, 1:n1)
        end
    end
    if n2 > 0
        if normalized
            I!(gram, rr)
        else
            gram[rr, rr] .= view(g.RAR, 1:n2, 1:n2)
        end
        gram[xr, rr] .= view(g.XAR, 1:n1, 1:n2)
        transpose!(view(gram, rr, xr), view(g.XAR, 1:n1, 1:n2))
    end
    if n3 > 0
        if normalized
            I!(gram, pr)
        else
            gram[pr, pr] .= view(g.PAP, 1:n3, 1:n3)
        end
        gram[rr, pr] .= view(g.RAP, 1:n2, 1:n3)
        gram[xr, pr] .= view(g.XAP, 1:n1, 1:n3)
        transpose!(view(gram, pr, rr), view(g.RAP, 1:n2, 1:n3))
        transpose!(view(gram, pr, xr), view(g.XAP, 1:n1, 1:n3))
    end
    return 
end

abstract type AbstractOrtho end
struct CholQR{TA} <: AbstractOrtho
    gramVBV::TA # to be used in view
end

function A_rdiv_B!(A, B::UpperTriangular)
    s = size(A, 2)
    A[:,1] .= view(A, :, 1) ./ B[1,1]
    for i in 2:s
        for j in 1:i-1
            A[:,i] .= view(A, :, i) - view(A, :, j) .* B[j,i]
        end
        A[:,i] .= view(A, :, i) ./ B[i,i]
    end
    return
end

realdiag!(M) = nothing
function realdiag!(M::AbstractMatrix{TC}) where TC <: Complex
    T = real(eltype(M))
    for i in 1:size(M, 1)
        M[i,i] = Complex(real(M[i,i]), zero(T))
    end
    return
end

function (ortho!::CholQR)(XBlocks::Blocks{Generalized}, sizeX = -1; update_AX=false, update_BX=false) where Generalized
    useview = sizeX != -1
    if sizeX == -1
        sizeX = size(XBlocks.block, 2)
    end
    X = XBlocks.block
    BX = XBlocks.B_block # Assumes it is premultiplied
    AX = XBlocks.A_block
    gram_view = view(ortho!.gramVBV, 1:sizeX, 1:sizeX)
    if useview
        Ac_mul_B!(gram_view, view(X, :, 1:sizeX), view(BX, :, 1:sizeX))
    else
        Ac_mul_B!(gram_view, X, BX)
    end
    realdiag!(gram_view)
    cholf = cholfact!(Hermitian(gram_view))
    R = cholf.factors
    if useview
        A_rdiv_B!(view(X, :, 1:sizeX), UpperTriangular(R))
        update_AX && A_rdiv_B!(view(AX, :, 1:sizeX), UpperTriangular(R))
        Generalized && update_BX && A_rdiv_B!(view(BX, :, 1:sizeX), UpperTriangular(R))
    else
        A_rdiv_B!(X, UpperTriangular(R))
        update_AX && A_rdiv_B!(AX, UpperTriangular(R))
        Generalized && update_BX && A_rdiv_B!(BX, UpperTriangular(R))
    end

    return 
end

struct LOBPCGIterator{Generalized, T, TA, TB, TL<:AbstractVector{T}, TR<:AbstractVector, TPerm<:AbstractVector{Int}, TV<:AbstractArray{T}, TBlocks<:Blocks{Generalized, T}, TO<:AbstractOrtho, TP, TC, TG, TM, TTrace}
    A::TA
    B::TB
    ritz_values::TL
    λperm::TPerm
    λ::TL
    V::TV
    residuals::TR
    largest::Bool
    XBlocks::TBlocks
    tempXBlocks::TBlocks
    PBlocks::TBlocks
    activePBlocks::TBlocks # to be used in view
    RBlocks::TBlocks
    activeRBlocks::TBlocks # to be used in view
    iteration::Base.RefValue{Int}
    currentBlockSize::Base.RefValue{Int}
    ortho!::TO
    precond!::TP
    constr!::TC
    gramABlock::TG
    gramBBlock::TG
    gramA::TV
    gramB::TV
    activeMask::TM
    trace::TTrace
end

"""
    LOBPCGIterator(A, X, largest::Bool, P=nothing, C=nothing) -> iterator

# Arguments

- `A`: linear operator;
- `X`: initial guess of the Ritz vectors; to be overwritten by the eigenvectors;
- `largest`: `true` if largest eigenvalues are desired and false if smallest;
- `P`: preconditioner of residual vectors, must overload `A_ldiv_B!`;
- `C`: constraint to deflate the residual and solution vectors orthogonal
    to a subspace; must overload `A_mul_B!`.
"""
LOBPCGIterator(A, X, largest::Bool, P=nothing, C=nothing) = LOBPCGIterator(A, nothing, X, largest, P, C)

"""
    LOBPCGIterator(A, B, X, largest::Bool, P=nothing, C=nothing) -> iterator

# Arguments

- `A`: linear operator;
- `B`: linear operator;
- `X`: initial guess of the Ritz vectors; to be overwritten by the eigenvectors;
- `largest`: `true` if largest eigenvalues are desired and false if smallest;
- `P`: preconditioner of residual vectors, must overload `A_ldiv_B!`;
- `C`: constraint to deflate the residual and solution vectors orthogonal
    to a subspace; must overload `A_mul_B!`;

"""
function LOBPCGIterator(A, B, X, largest::Bool, P=nothing, C=nothing)
    T = eltype(X)
    constr! = Constraint(C, B, X)
    precond! = RPreconditioner(P, X)

    nev = size(X, 2)
    if B isa Void
        XBlocks = Blocks(X, similar(X))
        tempXBlocks = Blocks(copy(X), similar(X))
        RBlocks = Blocks(similar(X), similar(X))
        activeRBlocks = Blocks(similar(X), similar(X))
        PBlocks = Blocks(similar(X), similar(X))
        activePBlocks = Blocks(similar(X), similar(X))
    else
        XBlocks = Blocks(X, similar(X), similar(X))
        tempXBlocks = Blocks(copy(X), similar(X), similar(X))
        RBlocks = Blocks(similar(X), similar(X), similar(X))
        activeRBlocks = Blocks(similar(X), similar(X), similar(X))
        PBlocks = Blocks(similar(X), similar(X), similar(X))
        activePBlocks = Blocks(similar(X), similar(X), similar(X))
    end
    ritz_values = zeros(T, nev*3)
    λ = zeros(T, nev)
    λperm = zeros(Int, nev*3)
    V = zeros(T, nev*3, nev*3)
    residuals = fill(real(T)(NaN), nev)
    iteration = Ref(1)
    currentBlockSize = Ref(nev)
    generalized = !(B isa Void)
    ortho! = CholQR(zeros(T, nev, nev))

    gramABlock = BlockGram(XBlocks)
    gramBBlock = BlockGram(XBlocks)

    gramA = zeros(T, 3*nev, 3*nev)
    gramB = zeros(T, 3*nev, 3*nev)

    activeMask = ones(Bool, nev)
    trace = LOBPCGTrace{Vector{real(T)},Vector{T}}()

    return LOBPCGIterator{generalized, T, typeof(A), typeof(B), typeof(λ), typeof(residuals), typeof(λperm), typeof(V), typeof(XBlocks), typeof(ortho!), typeof(precond!), typeof(constr!), typeof(gramABlock), typeof(activeMask), typeof(trace)}(A, B, ritz_values, λperm, λ, V, residuals, largest, XBlocks, tempXBlocks, PBlocks, activePBlocks, RBlocks, activeRBlocks, iteration, currentBlockSize, ortho!, precond!, constr!, gramABlock, gramBBlock, gramA, gramB, activeMask, trace)
end

function ortho_AB_mul_X!(blocks::Blocks, ortho!, A, B, bs=-1)
    # Finds BX
    bs == -1 ? B_mul_X!(blocks, B) : B_mul_X!(blocks, B, bs)
    # Orthonormalizes X and updates BX
    bs == -1 ? ortho!(blocks, update_BX=true) : ortho!(blocks, bs, update_BX=true)
    # Updates AX
    bs == -1 ? A_mul_X!(blocks, A) : A_mul_X!(blocks, A, bs)
    return 
end
function residuals!(iterator)
    sizeX = size(iterator.XBlocks.block, 2)
    A_mul_B!(iterator.RBlocks.block, iterator.XBlocks.B_block, Diagonal(view(iterator.ritz_values, 1:sizeX)))
    iterator.RBlocks.block .= iterator.XBlocks.A_block .- iterator.RBlocks.block
    # Finds residual norms
    for j in 1:size(iterator.RBlocks.block, 2)
        iterator.residuals[j] = 0
        for i in 1:size(iterator.RBlocks.block, 1)
            x = iterator.RBlocks.block[i,j]
            iterator.residuals[j] += real(x*conj(x))
        end
        iterator.residuals[j] = sqrt(iterator.residuals[j])
    end
    return
end

function update_mask!(iterator, residualTolerance)
    sizeX = size(iterator.XBlocks.block, 2)
    # Update active vectors mask
    iterator.activeMask .= view(iterator.residuals, 1:sizeX) .> residualTolerance
    iterator.currentBlockSize[] = sum(iterator.activeMask)
    return 
end

function update_active!(mask, bs::Int, blockPairs...)
    for (activeblock, block) in blockPairs
        activeblock[:, 1:bs] .= view(block, :, mask)
    end
    return
end

function precond_constr!(block, bs, precond!, constr!)
    precond!(view(block, :, 1:bs))
    # Constrain the active residual vectors to be B-orthogonal to Y
    constr!(view(block, :, 1:bs))
    return 
end
function block_grams_1x1!(iterator)
    # Finds gram matrix X'AX
    XAX!(iterator.gramABlock, iterator.XBlocks)
    return
end
function block_grams_2x2!(iterator, bs)
    sizeX = size(iterator.XBlocks.block, 2)
    #XAX!(iterator.gramABlock, iterator.XBlocks)
    XAR!(iterator.gramABlock, iterator.XBlocks, iterator.activeRBlocks, bs)
    RAR!(iterator.gramABlock, iterator.activeRBlocks, bs)
    XBR!(iterator.gramBBlock, iterator.XBlocks, iterator.activeRBlocks, bs)        
    iterator.gramABlock(iterator.gramA, view(iterator.ritz_values, 1:sizeX), sizeX, bs, 0)
    iterator.gramBBlock(iterator.gramB, sizeX, bs, 0, true)

    return
end
function block_grams_3x3!(iterator, bs)
    # Find R'AR, P'AP, X'AR, X'AP and R'AP
    sizeX = size(iterator.XBlocks.block, 2)
    #XAX!(iterator.gramABlock, iterator.XBlocks)
    XAR!(iterator.gramABlock, iterator.XBlocks, iterator.activeRBlocks, bs)
    XAP!(iterator.gramABlock, iterator.XBlocks, iterator.activePBlocks, bs)
    RAR!(iterator.gramABlock, iterator.activeRBlocks, bs)
    RAP!(iterator.gramABlock, iterator.activeRBlocks, iterator.activePBlocks, bs)
    PAP!(iterator.gramABlock, iterator.activePBlocks, bs)
    # Find X'BR, X'BP and P'BR
    XBR!(iterator.gramBBlock, iterator.XBlocks, iterator.activeRBlocks, bs)
    XBP!(iterator.gramBBlock, iterator.XBlocks, iterator.activePBlocks, bs)
    RBP!(iterator.gramBBlock, iterator.activeRBlocks, iterator.activePBlocks, bs)    
    # Update the gram matrix [X R P]' A [X R P]
    iterator.gramABlock(iterator.gramA, view(iterator.ritz_values, 1:sizeX), sizeX, bs, bs)
    # Update the gram matrix [X R P]' B [X R P]
    iterator.gramBBlock(iterator.gramB, sizeX, bs, bs, true)

    return
end

function sub_problem!(iterator, sizeX, bs1, bs2)
    subdim = sizeX+bs1+bs2
    if bs1 == 0
        gramAview = view(iterator.gramABlock.XAX, 1:subdim, 1:subdim)
        # Source of type instability
        realdiag!(gramAview)
        eigf = eigfact!(Hermitian(gramAview))
    else
        gramAview = view(iterator.gramA, 1:subdim, 1:subdim)
        gramBview = view(iterator.gramB, 1:subdim, 1:subdim)
        # Source of type instability
        realdiag!(gramAview)
        realdiag!(gramBview)
        eigf = eigfact!(Hermitian(gramAview), Hermitian(gramBview))
    end
    # Selects extremal eigenvalues and corresponding vectors
    selectperm!(view(iterator.λperm, 1:subdim), eigf.values, 1:subdim, rev=iterator.largest)
    iterator.ritz_values[1:sizeX] .= view(eigf.values, view(iterator.λperm, 1:sizeX))
    iterator.V[1:subdim, 1:sizeX] .= view(eigf.vectors, :, view(iterator.λperm, 1:sizeX))

    return
end

function update_X_P!(iterator::LOBPCGIterator{Generalized}, bs1, bs2) where Generalized
    sizeX = size(iterator.XBlocks.block, 2)
    x_eigview = view(iterator.V, 1:sizeX, 1:sizeX)
    r_eigview = view(iterator.V, sizeX+1:sizeX+bs1, 1:sizeX)
    p_eigview = view(iterator.V, sizeX+bs1+1:sizeX+bs1+bs2, 1:sizeX)
    r_blockview = view(iterator.activeRBlocks.block, :, 1:bs1)
    ra_blockview = view(iterator.activeRBlocks.A_block, :, 1:bs1)
    p_blockview = view(iterator.activePBlocks.block, :, 1:bs2)
    pa_blockview = view(iterator.activePBlocks.A_block, :, 1:bs2)
    if Generalized
        rb_blockview = view(iterator.activeRBlocks.B_block, :, 1:bs1)
        pb_blockview = view(iterator.activePBlocks.B_block, :, 1:bs2)
    end
    if bs1 > 0
        A_mul_B!(iterator.PBlocks.block, r_blockview, r_eigview)
        A_mul_B!(iterator.PBlocks.A_block, ra_blockview, r_eigview)
        if Generalized
            A_mul_B!(iterator.PBlocks.B_block, rb_blockview, r_eigview)
        end
    end
    if bs2 > 0
        A_mul_B!(iterator.tempXBlocks.block, p_blockview, p_eigview)
        A_mul_B!(iterator.tempXBlocks.A_block, pa_blockview, p_eigview)
        if Generalized
            A_mul_B!(iterator.tempXBlocks.B_block, pb_blockview, p_eigview)
        end
        iterator.PBlocks.block .= iterator.PBlocks.block .+ iterator.tempXBlocks.block
        iterator.PBlocks.A_block .= iterator.PBlocks.A_block .+ iterator.tempXBlocks.A_block
        if Generalized
            iterator.PBlocks.B_block .= iterator.PBlocks.B_block .+ iterator.tempXBlocks.B_block
        end
    end
    block = iterator.XBlocks.block
    tempblock = iterator.tempXBlocks.block
    A_mul_B!(tempblock, block, x_eigview)
    block = iterator.XBlocks.A_block
    tempblock = iterator.tempXBlocks.A_block
    A_mul_B!(tempblock, block, x_eigview)
    if Generalized
        block = iterator.XBlocks.B_block
        tempblock = iterator.tempXBlocks.B_block
        A_mul_B!(tempblock, block, x_eigview)
    end
    if bs1 > 0
        iterator.XBlocks.block .= iterator.tempXBlocks.block .+ iterator.PBlocks.block
        iterator.XBlocks.A_block .= iterator.tempXBlocks.A_block .+ iterator.PBlocks.A_block
        if Generalized
            iterator.XBlocks.B_block .= iterator.tempXBlocks.B_block .+ iterator.PBlocks.B_block
        end
    else
        iterator.XBlocks.block .= iterator.tempXBlocks.block
        iterator.XBlocks.A_block .= iterator.tempXBlocks.A_block
        if Generalized
            iterator.XBlocks.B_block .= iterator.tempXBlocks.B_block
        end
    end    
    return
end

function (iterator::LOBPCGIterator{Generalized})(residualTolerance, log) where {Generalized}
    sizeX = size(iterator.XBlocks.block, 2)
    iteration = iterator.iteration[]
    if iteration == 1
        iterator.constr!(iterator.XBlocks.block)
        ortho_AB_mul_X!(iterator.XBlocks, iterator.ortho!, iterator.A, iterator.B)
        # Finds gram matrix X'AX
        block_grams_1x1!(iterator)
        sub_problem!(iterator, sizeX, 0, 0)
        # Updates Ritz vectors X and updates AX and BX accordingly
        update_X_P!(iterator, 0, 0)
        residuals!(iterator)
        update_mask!(iterator, residualTolerance)
    elseif iteration == 2
        bs = iterator.currentBlockSize[]
        # Update active R blocks
        update_active!(iterator.activeMask, bs, (iterator.activeRBlocks.block, iterator.RBlocks.block))
        # Precondition and constrain the active residual vectors
        precond_constr!(iterator.activeRBlocks.block, bs, iterator.precond!, iterator.constr!)
        # Orthonormalizes R[:,1:bs] and finds AR[:,1:bs] and BR[:,1:bs]
        ortho_AB_mul_X!(iterator.activeRBlocks, iterator.ortho!, iterator.A, iterator.B, bs)
        # Find [X R] A [X R] and [X R]' B [X R]
        block_grams_2x2!(iterator, bs)
        # Solve the Rayleigh-Ritz sub-problem
        sub_problem!(iterator, sizeX, bs, 0)
        update_X_P!(iterator, bs, 0)
        residuals!(iterator)
        update_mask!(iterator, residualTolerance)
    else
        # Update active blocks
        bs = iterator.currentBlockSize[]
        # Update active R and P blocks
        update_active!(iterator.activeMask, bs, 
            (iterator.activeRBlocks.block, iterator.RBlocks.block), 
            (iterator.activePBlocks.block, iterator.PBlocks.block),
            (iterator.activePBlocks.A_block, iterator.PBlocks.A_block),
            (iterator.activePBlocks.B_block, iterator.PBlocks.B_block))
        # Precondition and constrain the active residual vectors
        precond_constr!(iterator.activeRBlocks.block, bs, iterator.precond!, iterator.constr!)
        # Orthonormalizes R[:,1:bs] and finds AR[:,1:bs] and BR[:,1:bs]
        ortho_AB_mul_X!(iterator.activeRBlocks, iterator.ortho!, iterator.A, iterator.B, bs)
        # Orthonormalizes P and updates AP
        iterator.ortho!(iterator.activePBlocks, bs, update_AX=true, update_BX=true)
        # Update the gram matrix [X R P]' A [X R P] and [X R P]' B [X R P]
        block_grams_3x3!(iterator, bs)
        # Solve the Rayleigh-Ritz sub-problem
        sub_problem!(iterator, sizeX, bs, bs)
        # Updates Ritz vectors X and updates AX and BX accordingly
        # And updates P, AP and BP
        update_X_P!(iterator, bs, bs)
        residuals!(iterator)
        update_mask!(iterator, residualTolerance)
    end
    if log
        return LOBPCGState(iteration, iterator.residuals[1:sizeX], iterator.ritz_values[1:sizeX])
    else
        return LOBPCGState(iteration, nothing, nothing)
    end
end

"""
The Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG)

Finds the k extremal eigenvalues and their corresponding eigenvectors satisfying `AX = λBX` returning `λ` and `X`, using `X` as an initial guess, where `k = size(X,2)`. The input `X` is overwritten with the solution.

`A` and `B` may be generic types but `Base.A_mul_B!(C, AorB, X)` and `AorB * X` must be defined for vectors and strided matrices `X` and `C`.

    lobpcg(A, [B,] largest, X0; kwargs...) -> results

# Arguments

- `A`: linear operator;
- `B`: linear operator;
- `largest`: `true` if largest eigenvalues are desired and false if smallest;
- `X0`: Initial guess, will not be modified. The number of columns is the number of eigenvectors desired.

## Keywords

- `log::Bool`: default is `false`; if `true`, `results.trace` will store iterations    
    states; if `false` only `results.trace` will be empty;

- `P`: preconditioner of residual vectors, must overload `A_ldiv_B!`;

- `C`: constraint to deflate the residual and solution vectors orthogonal
    to a subspace; must overload `A_mul_B!`;

- `maxiter`: maximum number of iterations; default is 200;

- `tol::Number`: tolerance to which residual vector norms must be under.

# Output

- `results`: a `LOBPCGResults` struct. `r.λ` and `r.X` store the eigenvalues and eigenvectors.

"""
function lobpcg(A, largest::Bool, X0::Union{AbstractMatrix, AbstractVector}; kwargs...)
    lobpcg(A, nothing, largest, X0; kwargs...)
end
function lobpcg(A, B, largest, X0;
                log=false, P=nothing, C=nothing, 
                tol=nothing, maxiter=200)

    X = copy(X0)
    T = eltype(X)
    n = size(X, 1)
    sizeX = size(X, 2)
    sizeX > n && throw("X column dimension exceeds the row dimension")

    iterator = LOBPCGIterator(A, B, X, largest, P, C)
    
    return lobpcg!(iterator, log=log, tol=tol, maxiter=maxiter)
end

"""
    lobpcg!(iterator::LOBPCGIterator; kwargs...) -> results
    
# Arguments

- `iterator::LOBPCGIterator`: a struct having all the variables required
    for the LOBPCG algorithm.

## Keywords

- `log::Bool`: default is `false`; if `true`, `results.trace` will store iterations    
    states; if `false` only `results.trace` will be empty;

- `maxiter`: maximum number of iterations; default is 200;

- `tol::Number`: tolerance to which residual vector norms must be under.

# Output

- `results`: a `LOBPCGResults` struct. `r.λ` and `r.X` store the eigenvalues and eigenvectors.

"""
function lobpcg!(iterator::LOBPCGIterator; log=false, tol=nothing, maxiter=200)
    T = eltype(iterator.XBlocks.block)
    X = iterator.XBlocks.block
    for j in 1:size(X,2)
        if all(x -> x==0, view(X, :, j))
            X[:,j] .= rand.()
        end
    end
    n = size(X, 1)
    sizeX = size(X, 2)
    residualTolerance = (tol isa Void) ? sqrt(eps(real(T))) : tol
    iterator.iteration[] = 1
    while iterator.iteration[] <= maxiter 
        state = iterator(residualTolerance, log)
        if log
            push!(iterator.trace, state)
        end
        iterator.currentBlockSize[] == 0 && break
        iterator.iteration[] += 1
    end
    iterator.λ .= view(iterator.ritz_values, 1:sizeX)

    results = LOBPCGResults(iterator.λ, X, iterator.iteration[], iterator.residuals[1:sizeX], residualTolerance, all((x)->(norm(x)<=residualTolerance), view(iterator.residuals, 1:sizeX)), maxiter, iterator.trace)

    return results
end
