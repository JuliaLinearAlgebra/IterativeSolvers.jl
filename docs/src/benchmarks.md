# Benchmarks


This document provides benchmarks of the iterative methods over well-known matrices. The matrices are all of size 10. For more benchmark
    options go to the [benchmark](https://github.com/JuliaLang/IterativeSolvers.jl/tree/master/benchmark)
    folder in the package.

```@contents
Pages = ["benchmarks.md"]
```

## Dense

### Kahan
```
                               Kahan (μs)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪ 4.842                                  
         jacobi │▪▪▪ 10.81                                
            sor │▪▪▪▪▪▪▪▪ 32.509                          
      chebyshev │▪▪▪ 10.021                               
           lsqr │▪▪ 6.944                                 
          gmres │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 119.711  
                └────────────────────────────────────────┘ 
```

### Vand
```
                                Vand (μs)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪▪ 46.338                                
         jacobi │▪▪▪▪▪▪▪▪ 156.701                         
            sor │▪▪▪ 58.549                               
      chebyshev │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 605.185  
           lsqr │ 7.835                                   
          gmres │▪▪▪▪▪▪ 121.579                           
                └────────────────────────────────────────┘ 
```

### Hilb
```
                                Hilb (μs)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪▪ 46.301                                
           ssor │ 7.61                                    
         jacobi │▪▪▪▪▪▪ 137.001                           
            sor │▪▪ 51.3                                  
             cg │ 11.302                                  
      chebyshev │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 736.975  
           lsqr │▪ 13.018                                 
          gmres │▪▪▪▪▪▪▪ 162.988                          
                └────────────────────────────────────────┘ 
```

### KMS
```
                                KMS (μs)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪▪▪▪ 14.64                               
           ssor │▪▪ 8.194                                 
         jacobi │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 120.137  
            sor │▪▪▪▪▪▪▪▪▪▪▪▪▪ 48.691                     
             cg │▪▪ 7.24                                  
      chebyshev │▪▪▪▪▪▪▪▪▪ 32.759                         
           lsqr │▪▪ 7.663                                 
          gmres │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 122.21  
                └────────────────────────────────────────┘ 
```

## Sparse

### Poisson
```
                              Poisson (μs)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪▪▪▪▪▪▪▪▪▪ 224.058                       
           ssor │▪▪▪▪ 104.173                             
         jacobi │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 567.76          
            sor │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 729.004  
             cg │ 7.524                                   
      chebyshev │▪▪▪ 76.862                               
           lsqr │▪ 12.84                                  
          gmres │▪▪▪▪▪ 122.64                             
                └────────────────────────────────────────┘ 
```

## Function

### SOLtest
```
                       SOLtest (μs)
         ┌────────────────────────────────────────┐ 
    lsqr │▪▪▪▪▪▪ 27.49                             
   gmres │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 134.085  
         └────────────────────────────────────────┘ 
```

