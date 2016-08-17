module Benchmarks

export addEqMethod, addEquation, buildSuite, runSuite

include(joinpath(Pkg.dir("IterativeSolvers"), "benchmark", "BenchmarkUtils.jl"))
using .BenchmarkUtils

using MatrixDepot
using BenchmarkTools
using IterativeSolvers


# Methods Tags
#
#   accessible        ->    Method accesses the linear operator's fields
#   inverse           ->    'A' Inverse must exist
#   symmetric         ->    'A' must be symmetric
#   pos-def           ->    'A' must be definite
#

# Linear Operators Tags
#
#   All methods tags  ->    'A' complies with method tag
#   MatrixDepot tags  ->    http://matrixdepotjl.readthedocs.io/en/latest/properties.html

# Set number of samples
BenchmarkTools.DEFAULT_PARAMETERS.samples = 400
# Set number of evaluations per sample to 1
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1

# Automatically generates an expression to call function.
# This only works for linear equation methods.
#
#   input:
#       method    ->  Linear equation method
#       traits    ->  Method's tags
#   output:
#       mcall     ->  Expression call to method
#       extra     ->  Generator of mcall variables
function buildCall(method::Function, traits)
    mcall = :()
    extra = Expr(:block)
    if method in [sor, ssor]
        mcall = :($method(A,b,ω))
        extra = quote
            ω = 0.5
        end
    elseif method == chebyshev  #chebyshev is a particular case
        mcall = :($method(A,b,λmin,λmax))
        extra = quote
            v = eigvals(full(A));
            mxv = maximum(v);
            mnv = minimum(v);
            λmin = mxv+(mxv-mnv)/100;
            λmax = mnv-(mxv-mnv)/100
        end
    else
        mcall = :($method(A,b))
    end
    return mcall, extra
end

# Defines an expression of a linear equation in a dictionary.
#
# mcall     ->  LE method call
# traits    ->  Method tags
# extra     ->  Generator mcall of variables
defEqMethod(mcall::Expr, traits, extra::Expr=Expr(:block)) =
Dict("mcall" => mcall, "extra" => extra, "traits" => Set(traits))

# Adds a linear equation definition to a dictionary
#
# dic           ->  Dictionary
# key           ->  Key
# definition    ->  Linear equation definition
# mcall         ->  LE method call
# traits        ->  Method tags
# extra         ->  Generator mcall of variables
# method    ->  Function that generates a matrix
addEqMethod(dic, key, definition::Dict) = (dic[key] = definition)
addEqMethod(dic, key, mcall::Expr, traits, extra::Expr=Expr(:block)) =
    addEqMethod(dic, key, defEqMethod(mcall,traits,extra))
function addEqMethod(dic, key, method::Function, traits)
  mcall, extra = buildCall(method, traits)
  addEqMethod(dic, key, mcall, traits, extra)
end

# Makes an equation definition
#
# path      ->  Path of the matrix on 'suite'
# traits    ->  Array of linear operator tags
# genA      ->  Matrix generator
# genb      ->  Vector generator
# dim       ->  Length of b, dim of A
function defEquation(path,traits,genA::Expr,genb::Expr= :(rand(size(A,2))) )
  A = Expr(:(=), :A, genA)
  b = Expr(:(=), :b, genb)
  Dict("path" => path, "traits" => Set(traits), "A" => A, "b" => b)
end
defEquation(path,traits,genA::Expr,dim::Integer; genb= :(rand($dim))) =
  defEquation(path,traits,genA,genb)

# Adds an equation definition to a dictionary
#
#   dic       ->  Dictionary
#   key       ->  Dictionary key
#   path      ->  Path of the matrix on 'suite'
#   traits    ->  Array of linear operator tags
#   genA      ->  Matrix generator
#   genb      ->  Vector generator
#   dim       ->  Length of b, dim of A
addEquation(dic, key, path, traits,genA::Expr,genb::Expr= :(rand(size(A,2)))) =
    dic[key] = defEquation(path,traits,genA,genb)
addEquation(dic, key, path, traits, genA::Expr, dim::Integer; genb= :(rand($dim))) =
  addEquation(dic,key,path,traits,genA,genb)

# Creates a path in a benchmarkgroup and returns the
# leaf of the new path.
function goToPath(bg::BenchmarkGroup, path)
    for key in path
        !(key in keys(bg)) && (bg[key] = BenchmarkGroup())
        bg = bg[key]
    end
    bg
end

# Wraps expression in a block, unless it already is.
function encapsule(expr::Expr)
  if expr.head != :block
    res = Expr(:block)
    push!(res.args, expr)
    return res
  end
  expr
end

# Builds the benchmarkable
function addCustomTest(suite,name,method,equation)
    path = equation["path"]
    traitsA = equation["traits"]
    genA = encapsule(equation["A"])
    genb = encapsule(equation["b"])
    bg = goToPath(suite, path)
    mcall = method["mcall"]
    traitsm = method["traits"]
    extra = method["extra"]
    !issubset(traitsm,traitsA) && throw(error("Method can't be used with equation"))
    init = Expr(:block)
    init.args = vcat(genA.args, genb.args, extra.args)
    bg[name] = eval(parse("@benchmarkable $mcall evals=1 samples=50 setup=($init)"))
end

# Links all the equation definitions with all the method definitions.
# Returns a BenchmarkGroup with all the test resulting from the cross product,
# ignoring imcompatibles.
#
#   methods     ->  Method definitions
#   equations   ->  Equation definitions
function buildSuite(methods, equations)
    suite = BenchmarkGroup()
    for equation in values(equations)
        for mkey in keys(methods)
            try
                addCustomTest(suite,mkey,methods[mkey],equation)
            end
        end
    end
    suite
end


###Methods Preconditions###
methods = Dict()

addEqMethod(methods, "jacobi", jacobi, ["inverse","accessible"])
addEqMethod(methods, "gauss_seidel", gauss_seidel, ["inverse","accessible"])
addEqMethod(methods, "sor", sor, ["inverse","accessible"])
addEqMethod(methods, "ssor", ssor, ["inverse","accessible", "symmetric"])
addEqMethod(methods, "cg", cg, ["inverse", "symmetric", "pos-def"])
addEqMethod(methods, "gmres", gmres, ["inverse"])
addEqMethod(methods, "lsqr", lsqr, ["inverse"])
addEqMethod(methods, "chebyshev", chebyshev, ["inverse", "accessible"])

# Custom call example
# addEqMethod(methods, "name", name, ["inverse", "accessible"], :(var1 = size(A); var2 = 200))

###Equations definitions###
equations = Dict()

#Dense matrix equations
addEquation(
    equations, "Hilb", ["Dense", "Hilb"],
    ["dense","inverse", "ill-cond", "symmetric", "pos-def", "accessible"],
    :(matrixdepot("hilb",10))
    )
addEquation(
    equations, "Kahan", ["Dense", "Kahan"],
    ["dense","inverse", "ill-cond", "accessible"],
    :(matrixdepot("kahan",10))
    )
addEquation(
    equations, "Vand", ["Dense", "Vand"],
    ["dense","inverse", "ill-cond", "accessible"],
    :(matrixdepot("vand",10))
    )
addEquation(
    equations, "KMS", ["Dense", "KMS"],
    ["dense","inverse", "ill-cond", "symmetric", "pos-def", "accessible"],
    :(matrixdepot("kms",10))
    )

#Sparse matrix equations
addEquation(
    equations, "Poisson", ["Sparse", "Poisson"],
    ["sparse","inverse", "symmetric", "pos-def", "eigen", "accessible"],
    :(matrixdepot("poisson",4))
    )

#Function matrix equations
addEquation(
    equations, "SOLtest", ["Function", "SOLtest"],
    ["function","inverse"],
    :(buildSol(10)),
    10
    )

suite = buildSuite(methods, equations)

end
