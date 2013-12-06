export PosSemidefException

type PosSemidefException <: Exception
    msg :: UTF8String
    PosSemidefException(msg::String="Matrix was not positive semidefinite") = new(msg)
end

