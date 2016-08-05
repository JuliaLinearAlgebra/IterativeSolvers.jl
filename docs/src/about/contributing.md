# Contributing

Contributions are always welcome, as are feature requests and suggestions. Feel
free to open [issues](https://help.github.com/articles/creating-an-issue/) and [pull requests](https://help.github.com/articles/creating-a-pull-request/) at any time. If you aren't familiar with git or Github please start [now](https://help.github.com/articles/good-resources-for-learning-git-and-github/).

## Setting workspace up

Julia's internal package manager makes it easy to install and modify packages
from Github. Any package hosted on Github can be installed via `Pkg.clone` by
providing the [repository's URL](https://help.github.com/articles/which-remote-url-should-i-use/), so installing a fork on your system is a simple task.

```julia
Pkg.clone("https://github.com/johndoe/IterativeSolvers.jl")
```

It is to note here if you have the original package installed the fork will
replace it, this is not a problem.

Now find your fork's location

```julia
Pkg.dir("IterativeSolvers")
```

Once there you will notice you are on the master branch, whenever a package is
imported Julia will use the code in the current branch, this means checking out
other git branches will let you use/test whatever there is.
