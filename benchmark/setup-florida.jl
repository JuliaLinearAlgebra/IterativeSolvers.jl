#!/usr/bin/env julia

#Benchmark iterative methods against matrices from the University of Florida
#sparse collection

VERSION ≥ v"0.4" || error("Julia 0.4 required")

#Root URL to the matrix collection
UFL_URL_ROOT = "http://www.cise.ufl.edu/research/sparse/matrices"
#Plain text file containing list of matrices
MASTERLIST = "matrices.txt"
#Download UFL collection to this directory
BASEDIR = "florida"

# 1. Read in master list of matrices from BASEDIR/MASTERLIST
#    If absent, generate this file. Requires Gumbo.jl
Pkg.installed("Gumbo")==nothing || using Gumbo
isdir(BASEDIR) || mkdir(BASEDIR)
if !isfile(joinpath(BASEDIR, MASTERLIST))
    Pkg.installed("Gumbo")===nothing && error("Parsing list from UFL website requires Gumbo.jl")
    #Parse the list on the UFL website. Requires the Gumbo HTML parser

    dlfilename = joinpath(BASEDIR, "florida.html")
    isfile(dlfilename) || download(UFL_URL_ROOT*"/list_by_id.html", dlfilename)

    doc = parsehtml(open(readall, dlfilename))

    matrices = Tuple{AbstractString, AbstractString}[]
    for elem in preorder(doc.root)
        if isa(elem, HTMLElement) && tag(elem)==:a
            relurl = getattr(elem, "href")
            if endswith(relurl, ".mat")
                tokens = split(relurl, '/')
                group = tokens[end-1]
                name = split(tokens[end], '.')[1]

                push!(matrices, (group, name))
            end
        end
    end

    #Write master list
    open(joinpath(BASEDIR, MASTERLIST), "w") do io
        for (group, name) in matrices
            println(io, group, '\t', name)
        end
    end
else #Read list of matrices from saved file
    matrices = Tuple{AbstractString, AbstractString}[]
    open(joinpath(BASEDIR, MASTERLIST)) do io
        for l in eachline(io)
            tokens = split(l)
            push!(matrices, (tokens[1], tokens[2]))
        end
    end
end

info("Listing contains $(length(matrices)) matrices")

# 2. Download all the matrices
for (group, matrix) in matrices
    groupdir = joinpath(BASEDIR, group)
    isdir(groupdir) || mkdir(groupdir)

    if !isfile(joinpath(groupdir, matrix*".mat"))
        info("Downloading $group/$matrix.mat")
        download(joinpath(UFL_URL_ROOT, "..", "mat", group, matrix*".mat"),
            joinpath(BASEDIR, group, matrix*".mat"))
    end
end
info("Downloaded $(length(matrices)) matrices")


