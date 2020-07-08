# SimpleRoots

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/SimpleRoots.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/SimpleRoots.jl/dev)
[![Build Status](https://travis-ci.com/rafaqz/SimpleRoots.jl.svg?branch=master)](https://travis-ci.com/rafaqz/SimpleRoots.jl)
[![codecov.io](http://codecov.io/github/rafaqz/SimpleRoots.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/SimpleRoots.jl?branch=master)


A light-weight collection of simple root-finding algorithms.

These are intended to be simple but high-performance routines 
with no allocation or construction overheads. 
They also work with Unitful.jl, or any other types that define 
basic math and comparison operations. 

The package code and tests are small and easy to understand, 
and are guaranteed to stay that way.

Included are bracketed `findzero` methods including:
- Brent
- Bisection
- Secant

A tuple is returned, containing the value and a `Bool` for sucsess:

```julia-repl
julia> findzero(sin, Secant(-0.5, 0.5))
(0.0, true)
```

A basic quadratic solver `quad` is also included:

```julia-repl
julia> quad(1.0, 3.0, -4.0)
(-4.0, 1.0)
```

Accuracy and quality of methods is probably lower, and is far less
well tested than Roots.jl. For non-modelling purposes, use Roots.jl.
