module BenchmarkUtils

export getunit, prettytime, buildSol

using IterativeSolvers

#Functions

# Get closest unit to t
function getunit(t)
    t < 1e3 && (return "ns")
    t < 1e6 && (return "μs")
    t < 1e9 && (return "ms")
    "s"
end

# Gets time based on unit
function prettytime(t,unit)
    unit == "ns" && (return t)
    unit == "μs" && (return t/1e3)
    unit == "ms" && (return t/1e6)
    t/1e9
end

function buildSol(dim)
  fmul  = (out, b) -> Aprodxxx!(out,b,1,dim,dim)
  fcmul = (out, b) -> Aprodxxx!(out,b,2,dim,dim)
  MatrixCFcn{Int}(dim, dim, fmul, fcmul)
end

function Aprodxxx!(y, x, mode, m, n )
  # if mode = 1, computes y = A*x
  # if mode = 2, computes y = A'*x
  # for some matrix  A.
  #
  # This is a simple example for testing  LSQR.
  # It uses the leading m*n submatrix from
  # A = [ 1
  #       1 2
  #         2 3
  #           3 4
  #             ...
  #               n ]
  # suitably padded by zeros.
  #
  # 11 Apr 1996: First version for distribution with lsqr.m.
  #              Michael Saunders, Dept of EESOR, Stanford University.

  if mode == 1
      y[1] = x[1]
      for i = 2:n
          y[i] = i*x[i] + (i-1)*x[i-1]
      end
      for i = n+1:m
          y[i] = 0
      end
  else
      mn = min(m, n)
      for i = 1:mn-1
          y[i] = i*(x[i]+x[i+1])
      end
      y[mn] = mn*x[mn]
      for i = m+1:n
          y[i] = 0
      end
  end
  y
end

end
