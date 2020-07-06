module SimpleRoots

export Bisection, Secant, Brent 

export Lower, Upper, Both

export findzero, quad

abstract type AbstractRootFinder end

abstract type AbstractBisection <: AbstractRootFinder end

bracket(method::AbstractBisection) = method.bracket

struct Bisection{B} <: AbstractBisection 
    bracket::B
end
Bisection(arg1, arg2) = Bisection((arg1, arg2))

"""
    Secant(bracket)
    Secant(bracket...)

Secant method of bracketed root-finding.

Bracket can be a length 2 vector, a tuple or 2 argument
"""
struct Secant{B} <: AbstractBisection 
    bracket::B
end
Secant(arg1, arg2) = Secant((arg1, arg2))

struct Brent{B} <: AbstractBisection 
    bracket::B
end
Brent(arg1, arg2) = Brent((arg1, arg2))

"""
    findzero(f, bracket; atol=1e-7, max_iter=100)
    findzero(f, method::AbstractBisection; atol=1e-7, max_iter=100)

Find roots using the passed in method. If a bracket is passed in without
specifying the method, use the [`Secant`](@ref) method.
"""
findzero(f, bracket; kwargs...) = findzero(f, Secant(bracket); kwargs...)
findzero(f, method::AbstractRootFinder; atol=1e-7, max_iter=100) =
    findzero(f, method, atol, max_iter)

"""
    findzero(f, method::Secant, atol, max_iter)
    findzero(f, method::Secant; atol=1e-7, max_iter=100)

Find root using the secant method.
"""
function findzero(f, method::Secant, atol, max_iter)
    local x
    x0, x1 = bracket(method)
    y0, y1 = f(x0), f(x1) 

    for _ in 1:max_iter
        x = x1 - y1 * (x1 - x0) / (y1 - y0)
        if abs(x - x1) < atol 
            return x, true 
        end
        x0 = x1
        y0 = y1
        x1 = x
        y1 = f(x1)
    end
    return x, false
end

"""
    findzero(f, method::Bisection, atol, max_iter)
    findzero(f, method::Bisection; atol=1e-7, max_iter=100)

Find root using the bisection method.
"""
function findzero(f, method::Bisection, atol, max_iter)
    a, b = bracket(method)
    fa = f(a)
    local c
    for _ in 1:max_iter
        if (b - a) <= atol 
            return b, true 
        end

        c = (a + b)/2 
        fc = f(c)
        if fa * fc > zero(fa * fc)
            a = c # Root is in the right half of [a, b].
            fa = fc
        elseif fc == zero(fc)
            return c, true
        else
            b = c # Root is in the left half of [a, b].
        end
    end

    return b, false
end

"""
    findzero(f, method::Brent, atol, maxiter)
    findzero(f, method::Brent; atol=1e-7, max_iter=100)

Find root using Brents method. Returns a tuple of the found value
and a Bool specifying wether the root was found within the tolerance, or not.
"""
function findzero(f, method::Brent, atol, maxiter)
    local e, d
    epsval = eps()
    a, b = bracket(method)
    fa, fb = f(a), f(b)
    c = b
    fc = fb
    for iter = 1:maxiter
        if (fb > zero(fb) && fc > zero(fc)) || (fb < zero(fb) && fc < zero(fc))
            c = a
            fc = fa
            d = b - a
            e = d
        end
        if abs(fc) < abs(fb)
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end
        tol1 = 2epsval * abs(b) + atol/2
        xm = (c - b)/2
        if (abs(xm) <= tol1 || fb == 0fb)
            return b, true
        end
        if abs(e) >= tol1 && abs(fa) > abs(fb)
            s = fb / fa
            if (a == c)
                p = 2xm * s
                q = 1 - s
            else
                q = fa / fc
                r = fb / fc
                p = s * (2xm * q * (q - r) - (b - a) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)
            end
            if p > zero(p) 
                q = -q 
            end
            p = abs(p)
            if (2p < min(3xm * q - abs(tol1 * q), abs(e * q)))
                e = d
                d = p / q
            else
                d = xm
                e = d
            end
        else
            d = xm
            e = d
        end
        a = b
        fa = fb
        if (abs(d) > tol1)
            b = b + d
        else
            b = b + tol1 * sign(xm)
        end
        fb = f(b)
    end

    return b, false
end


"""
Used to specifiy the result of quad, either
[`Lower`](@ref), [`Upper`](@ref) or [`Both`](@ref)
"""
abstract type QuadraticResult end

"""
Return the lower result value
"""
struct Lower <: QuadraticResult end

"""
Return the upper result value
"""
struct Upper <: QuadraticResult end

"""
Return both lower and upper result values as a tuple
"""
struct Both <: QuadraticResult end

struct QuadraticError <: Exception 
   var::String
end

Base.showerror(io::IO, e::QuadraticError) = print(io, e.var, "!")

"""
    function quad(res::QuadraticResult, a, b, c)

Simple quadratic equation solver.

Returns the result specified with `res`
[`Lower`](@ref), [`Upper`](@ref) or [`Both`](@ref)
"""
function quad(res::QuadraticResult, a, b, c)
    x = (b^2 - 4 * a * c) 
    x < zero(x) && throw(QuadraticError("Imaginary roots in quadratic"))
    if a == zero(a)
        if b == zero(b)
            c != zero(c) && throw(QuadraticError("Can't solve quadratic"))
            # return zero with correct units
            return out(res, zero(-c / b))
        else
            return out(res, -c / b)
        end
    else 
        return side(res, x, a, b)
    end
end
quad(a, b, c) = quad(Both(), a, b, c)

@inline side(::Both, x, a, b) = side(Lower(), x, a, b), side(Upper(), x, a, b) 
@inline side(::Upper, x, a, b) = (-b + sqrt(x)) / 2a
@inline side(::Lower, x, a, b) = (-b - sqrt(x)) / 2a

@inline out(::Both, x) = (x, x)
@inline out(::Upper, x) = x
@inline out(::Lower, x) = x

end # module
