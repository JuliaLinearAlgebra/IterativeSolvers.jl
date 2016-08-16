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
IterativeSolvers.PlainHistory
IterativeSolvers.RestartedHistory
```

**`Functions`**

```@docs
IterativeSolvers.nextiter!
IterativeSolvers.reserve!
IterativeSolvers.shrink!
IterativeSolvers.setmvps
IterativeSolvers.setmtvps
IterativeSolvers.setconv
IterativeSolvers.showplot
```

## KrylovSubspace Internals

**`Functions`**

```@docs
IterativeSolvers.lastvec
IterativeSolvers.nextvec
IterativeSolvers.init!
IterativeSolvers.initrand!(::KrylovSubspace)
IterativeSolvers.appendunit!
IterativeSolvers.orthogonalize
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
IterativeSolvers.update!
IterativeSolvers.initrand!(::Vector)
```
