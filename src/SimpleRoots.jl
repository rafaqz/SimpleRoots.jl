__precompile__()

module SimpleRoots

export Bisection, Secant, ZBrent, FalsePosition, find_zero, quadm, quadp

abstract type AbstractBisection end
type Bisection <: AbstractBisection end
type Secant <: AbstractBisection end
type FalsePosition <: AbstractBisection end
type ZBrent <: AbstractBisection end

function find_zero(f, bracket; atol=1e-7, max_iter=100)
    find_zero(Secant(), f, bracket; atol, max_iter)
end

function find_zero(f, bracket, ::Secant, atol, max_iter)
    local x
    x0, x1 = bracket
    y0, y1 = f(x0), f(x1) 

    for _ in 1:max_iter
        x = x1 - y1 * (x1 - x0) / (y1 - y0)
        if abs(x - x1) < atol return x end
        x0 = x1
        y0 = y1
        x1 = x
        y1 = f(x1)
    end
    error("Root not found")
end

function find_zero(f, bracket, ::Bisection, atol, max_iter)
    a, b = bracket
    fa = f(a)
    local c
    for _ in 1:max_iter
        if (b - a) <= atol return b end

        c = (a + b)/2 
        fc = f(c)
        if fa * fc > 0(fa * fc)
            a = c # Root is in the right half of [a, b].
            fa = fc
        elseif fc == 0fc
            return c
        else
            b = c # Root is in the left half of [a, b].
        end
    end

    error("Root not found")
end

function find_zero(f, bracket, ::ZBrent; tolz=1e-7, maxiter=30)
    local e, d
    epsval = eps()
    a, b = bracket
    fa, fb = f(a), f(b)
    c = b
    fc = fb
    for iter = 1:maxiter
        if (fb > 0fb && fc > 0fc) || (fb < 0fb && fc < 0fc)
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
        tol1 = 2epsval * abs(b) + tolz/2
        xm = (c - b)/2
        if (abs(xm) <= tol1 || fb == 0fb)
            return b
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
            if p > 0p q = -q end
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

    return b
end


"""
Solves quadratic equations - finds smaller root.
"""
function quadm(a, b, c)
    x = (b * b - 4.0 * a * c) 
    if x < 0.0x
        error("imaginary roots in quadratic")
    end
    if a == 0.0a
        if b == 0.0b
            quadm = 0.0
            if (c != 0.0c) error("error: cant solve quadratic") end
        else
            quadm = -c / b
        end
    else
        quadm = (-b - sqrt(x)) / (2.0 * a)
    end

    return quadm
end

"""
Solves quadratic equations - finds larger root.
"""
function quadp(a, b, c)
    x = (b * b - 4.0 * a * c)
    if x < 0.0x
        error("imaginary roots in quadratic")
    end
    if a == 0.0a
        if b == 0.0b
            quadp = 0.0
            if (c != 0.0c) error("error: cant solve quadratic") end
        else
            quadp = -c / b
        end
    else
        quadp = (-b + sqrt(x)) / (2.0 * a)
    end
    return quadp
end

end # module
