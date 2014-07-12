#Test orthogonalization methods
using Base.Test

debug = false
srand(1)

m=n=3
A=randn(m, n) |> float32

if v"0.2" <= VERSION < v"0.3-"
    QRhouseholder = qr(A, false)
else
    QRhouseholder = qr(A, thin=false)
end
debug && @show QRhouseholder

for algorithm in (cgs, cmgs)
	for normalization in (norm_none, norm_naive, norm_pythag)
		for reorth_criterion in (never, always, rutishauser(√2), giraudlangou(10.0))
			debug && println(algorithm(nothing, normalization, ReorthogonalizationAlg(reorth_criterion)))
			QRgs = qrfact!(copy(A), algorithm(nothing, normalization, ReorthogonalizationAlg(reorth_criterion)))
			@test_approx_eq QRgs.Q*QRgs.R A
			normalization == norm_none && continue
			debug && @show QRhouseholder[1], QRgs.Q
			@test abs(norm(QRhouseholder[1] - QRgs.Q) - iround(norm(QRhouseholder[1] - QRgs.Q))) < √eps(eltype(A))
			debug && @show QRhouseholder[2], QRgs.R
			@test_approx_eq abs(QRhouseholder[2]) abs(QRgs.R)
		end
	end
end
