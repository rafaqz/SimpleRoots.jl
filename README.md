# SimpleRoots

A light-weight collection of simple root-finding algorithms.

These are intended to be simple but high-performance routines 
with no allocation or construction overheads. 
They also work with Unitful.jl any other types that define 
basic math and comparison operations. 

The package code and tests are small and easy to understand, 
and are guaranteed to stay that way.

Included are bracketed `findzero` methods including:
- Brent
- Bisection
- Secant


```julia-repl
julia> findzero(cos, Secant(0.0, pi))
(0.0, true)
```

And a basic quadratic solver `quad`:

```julia-repl
julia> quad(1.0, 3.0, -4.0)
(-4.0, 1.0)
```

Accuracy and quality of methods is definitely lower and fas less
well tested than Roots.jl. For non-modelling purposes, use Roots.jl.
