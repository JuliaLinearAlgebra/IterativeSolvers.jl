#Test orthogonalization methods
using Base.Test

debug = false
m=n=4
A=randn(m, n)

srand(1)

QRhouseholder = qr(A, thin=false)

for algorithm in (cgs, cmgs), normalization in (norm_none, norm_naive, norm_pythag)
	debug && println(algorithm(nothing, normalization))
	QRgs = qrfact!(copy(A), algorithm(nothing, normalization))
	@test_approx_eq QRgs.Q*QRgs.R A
	normalization == norm_none && continue
	debug && @show QRhouseholder[1], QRgs.Q
	debug && @show norm(QRhouseholder[1] - QRgs.Q) - iround(norm(QRhouseholder[1] - QRgs.Q)), m*n*eps()
	@test norm(QRhouseholder[1] - QRgs.Q) - iround(norm(QRhouseholder[1] - QRgs.Q)) < m*n*eps()
	debug && @show QRhouseholder[2], QRgs.R
	@test_approx_eq abs(QRhouseholder[2]) abs(QRgs.R)
end
