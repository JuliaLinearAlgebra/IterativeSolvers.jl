# [IDR(s)](@id IDRs)

The Induced Dimension Reduction method is a family of simple and fast Krylov subspace algorithms for solving large nonsymmetric linear systems. The idea behind the IDR(s) variant is to generate residuals that are in the nested subspaces of shrinking dimensions.

## Usage

```@docs
idrs
idrs!
```

## Implementation details

The current implementation is based on the [MATLAB version](http://ta.twi.tudelft.nl/nw/users/gijzen/idrs.m) by 
Van Gijzen and Sonneveld. For background see [^Sonneveld2008], [^VanGijzen2011], [the IDR(s) webpage](http://homepage.tudelft.nl/1w5b5/idrs-software.html)
and the IDR chapter in [^Meurant2020]. 

!!! tip
    IDR(s) can be used as an [iterator](@ref Iterators).

[^Sonneveld2008]: IDR(s): a family of simple and fast algorithms for solving large nonsymmetric linear systems. P. Sonneveld and M. B. van Gijzen SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035--1062, 2008
[^VanGijzen2011]: Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits Bi-orthogonality Properties. M. B. van Gijzen and P. Sonneveld ACM Trans. Math. Software, Vol. 38, No. 1, pp. 5:1-5:19, 2011
[^Meurant2020]: The IDR family. G. Meurant and J. Duintjer Tebbens. In: Krylov Methods for Nonsymmetric Linear Systems. Springer Series in Computational Mathematics, vol 57. Springer, 2020. [doi:10.1007/978-3-030-55251-0_10](https://doi.org/10.1007/978-3-030-55251-0_10)

