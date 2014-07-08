#Test orthogonalization methods
using Base.Test

debug = false
m=n=4
A=randn(m, n)

QRhouseholder = qr(A, thin=false)

for algorithm in (cgs, cmgs)
	debug && println(algorithm)
	QRgs = qrfact!(copy(A), algorithm)
	@test norm(QRhouseholder[1] - QRgs.Q) - iround(norm(QRhouseholder[1] - QRgs.Q)) < m*n*eps()
	@test_approx_eq abs(QRhouseholder[2]) abs(QRgs.R)
end
