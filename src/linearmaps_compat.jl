import LinearAlgebra: mul!
using LinearMaps

mul!(y, A::LinearMaps.WrappedMap, x) = A_mul_B!(y, A, x)
