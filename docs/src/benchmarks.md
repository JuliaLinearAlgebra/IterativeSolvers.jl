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
   gauss_seidel │▪ 4.871                                  
         jacobi │▪▪▪ 11.12                                
            sor │▪▪▪▪▪▪▪▪▪▪▪ 40.681                       
      chebyshev │▪▪ 8.954                                 
           lsqr │▪▪ 7.041                                 
          gmres │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 115.271  
                └────────────────────────────────────────┘ 
```

### Vand
```
                                Vand (μs)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪▪ 47.566                                
         jacobi │▪▪▪▪▪▪▪ 133.11                           
            sor │▪▪▪ 51.459                               
      chebyshev │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 597.899  
           lsqr │ 7.241                                   
          gmres │▪▪▪▪▪▪ 122.333                           
                └────────────────────────────────────────┘ 
```

### Hilb
```
                                Hilb (ms)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪▪ 0.048031                              
           ssor │ 0.007759                                
         jacobi │▪▪▪▪▪▪ 0.176405                          
            sor │▪▪ 0.055434                              
             cg │ 0.006535                                
      chebyshev │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 1.0148  
           lsqr │ 0.010512                                
          gmres │▪▪▪▪▪ 0.145835                           
                └────────────────────────────────────────┘ 
```

### KMS
```
                                KMS (μs)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪▪▪▪ 14.992                              
           ssor │▪▪ 7.684                                 
         jacobi │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 115.291      
            sor │▪▪▪▪▪▪▪▪▪▪▪▪ 49.651                      
             cg │▪▪ 7.565                                 
      chebyshev │▪▪▪▪▪▪▪▪ 32.077                          
           lsqr │▪▪ 10.237                                
          gmres │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 131.536  
                └────────────────────────────────────────┘ 
```

## Sparse

### Poisson
```
                              Poisson (μs)
                ┌────────────────────────────────────────┐ 
   gauss_seidel │▪▪▪▪▪▪▪▪▪▪ 229.247                       
           ssor │▪▪▪▪▪▪ 135.41                            
         jacobi │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 575.374         
            sor │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 743.994  
             cg │ 11.721                                  
      chebyshev │▪▪▪ 82.643                               
           lsqr │▪ 13.447                                 
          gmres │▪▪▪▪▪ 129.751                            
                └────────────────────────────────────────┘ 
```

## Function

### SOLtest
```
                       SOLtest (μs)
         ┌────────────────────────────────────────┐ 
    lsqr │▪▪▪▪▪▪▪▪ 32.106                          
   gmres │▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ 121.715  
         └────────────────────────────────────────┘ 
```

