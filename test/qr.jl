#Test orthogonalization methods
using Base.Test

m=n=4
A=randn(m, n)

QRhouseholder = qr(A, thin=false)

QRcgs = qrfact!(A, cgs)

@test norm(QRhouseholder[1] - QRcgs.Q) % 1 < m*n*eps()
@test_approx_eq abs(QRhouseholder[2]) abs(QRcgs.R)

