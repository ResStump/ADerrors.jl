# ADerrors.jl

Error propagation and analysis of Monte Carlo data with the (``\Gamma``) method and automatic differentiation in `Julia`

The full documentation of the package is available via the usual
[Julia `REPL` help
mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Help-mode-1) and
online in [HTML format](https://ific.uv.es/~alramos/docs/ADerrors/).

This work is an implementation of several ideas in data analysis. If you use this package for your scientific work, please consider citing:
- U. Wolff, "Monte Carlo errors with less errors".
  Comput.Phys.Commun. 156 (2004) 143-153. DOI: 10.1016/S0010-4655(03)00467-3
- F. Virotta, "Critical slowing down and error analysis of lattice QCD simulations." PhD thesis.
- Stefan Schaefer, Rainer Sommer, Francesco Virotta, "Critical slowing
  down and error analysis in lattice QCD simulations". Nucl.Phys.B 845 (2011) 93-119.
- A. Ramos, "Automatic differentiation for error analysis of Monte Carlo data". Comput.Phys.Commun. 238 (2019) 19-35. DOI: 10.1016/j.cpc.2018.12.020. 
- M. Bruno, R. Sommer, In preparation.

## Installation

The package in not in the general registry. Still one can use the package manager
```julia
julia> import Pkg
(v1.1) pkg> add https://gitlab.ift.uam-csic.es/alberto/aderrors.jl
```

## Tutorial

It is better to start with the [Getting started](https://ific.uv.es/~alramos/docs/ADerrors/tutorial/) guide.



