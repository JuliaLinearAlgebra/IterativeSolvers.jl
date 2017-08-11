# Internal Documentation

Documentation for `IterativeSolvers.jl`'s internals.

```@contents
Pages = ["internal.md"]
Depth = 4
```

## Index

```@index
Pages = ["internal.md"]
```

## ConvergenceHistory Internals

**`Typealiases`**

```@docs
IterativeSolvers.RestartedHistory
```

**`Functions`**

```@docs
IterativeSolvers.nextiter!
IterativeSolvers.reserve!
IterativeSolvers.shrink!
IterativeSolvers.setconv
IterativeSolvers.showplot
```

## Other Functions


```@docs
IterativeSolvers.idfact
IterativeSolvers.isconverged
IterativeSolvers.thickrestart!
IterativeSolvers.harmonicrestart!
IterativeSolvers.plot_collection
IterativeSolvers.plotable
IterativeSolvers.Adivtype
IterativeSolvers.Amultype
IterativeSolvers.randx
IterativeSolvers.zerox
IterativeSolvers.initrand!(::Vector)
```
