###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    ADerrorsUtils.jl
### created: Fri Jun 26 19:30:27 2020
###                               

"""
    root_error(fnew::Function, x0::Float64, data::Vector{uwreal}) 

Returns `x` such that `fnew(x, data) = 0` (i.e. a root of the function). The error in `data` is propagated to the root position `x`.
```@example
using ADerrors # hide
# First define some arbitrary data
data = Vector{uwreal}(undef, 3)
data[1] = uwreal([1.0, 0.2],   120)
data[2] = uwreal([1.2, 0.023], 121)
data[3] = uwreal(rand(1000),   122)

# Now define a function
f(x, p) = x + p[1]*x + cos(p[2]*x+p[3])

# Find its root using x0=1.0 as initial
# guess of the position of the root
x = root_error(f, 1.0, data)
uwerr(x)
println("Root: ", x)

# Check
z = f(x, data)
uwerr(z)
println("Better be zero (with zero error): ", z)
```
"""
function root_error(fnew::Function, x::Float64,
                    data::Vector{uwreal}) 

    function fvec(x::Vector)
        return fnew(x[1], x[2:end])
    end
    xv = zeros(Float64, length(data)+1)
    for i in 1:length(data)
        xv[i+1] = data[i].mean
    end

    xdt   = @view xv[2:end]
    f     = x -> fnew(x, xdt)
    D(f)  = x -> ForwardDiff.derivative(f,float(x))
    xv[1] = Roots.find_zero((f,D(f)), x, Roots.Newton()) 

    cfg   = GradientConfig(fvec, xv, Chunk{2}());
    grad  = ForwardDiff.gradient(fvec, xv, cfg)
    gw    = @view grad[2:end]

    return addobs(data, -gw ./ grad[1], xv[1])
end

function chiexp(hess::Array{Float64, 2}, data::Vector{uwreal}, W::Vector{Float64})

    m = length(data)
    n = size(hess, 1) - m

    hm = view(hess, 1:n, n+1:n+m)
    sm = Array{Float64, 2}(undef, n, m)
    for i in 1:n, j in 1:m
        sm[i,j] = hm[i,j] / sqrt.(W[j])
    end
    maux = sm * sm'
    hi   = LinearAlgebra.pinv(maux)
    Px   = -hm' * hi * hm

    for i in 1:m
        Px[i,i] = W[i] + Px[i,i]
    end

    return trcov(Px, data)
end
    
function chiexp(hess::Array{Float64, 2}, data::Vector{uwreal}, W::Array{Float64, 2})

    m = length(data)
    n = size(hess, 1) - m
    
    Lm = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(W))
    Li = LinearAlgebra.inv(Lm.L)
    
    hm = view(hess, 1:n, n+1:n+m)
    sm = hm * Li'
    
    maux = sm * sm'
    hi   = LinearAlgebra.pinv(maux)
    Px   = W - hm' * hi * hm

    return trcov(Px, data)
end

@doc raw"""
    chiexp(chisq::Function,
                xp::Vector{Float64}, 
                data::Vector{uwreal};
                W = Vector{Float64}())

Given a ``\chi^2(p, d)``, function of the fit parameters `p[:]` and the data `d[:]`, compute the expected value of the ``\chi^2(p, d)``.

### Arguments

- `chisq`: Must be a function of two vectors (i.e. `chisq(p::Vector, d::Vector)`). The function is assumed to have the form

``\chi^2(p, d) = \sum_i [d_i - f_i(p)]W_{ij}[d_j - f_j(p)]``

where the function ``f_i(p)`` is an arbitrary function of the fit parameters. In simple words, the function is assumed to be quadratic in the data.
- `xp`: A vector of `Float64`. The value of the fit parameters at the minima.
- `data`: A vector of `uwreal`. The data whose fluctuations enter in the evaluation of the `chisq`.
- `W`: A matrix. The weights that enter in the evaluation of the `chisq` function. If a vector is passed, the matrix is assumed to be diagonal (i.e. **uncorrelated** fit). If no weights are passed, the routines assumes that `W` is diagonal with entries given by the inverse errors squared of the data (i.w. the `chisq` is weighted with the errors of the data). 

### Example 
```@example
using ADerrors # hide

```

"""
function chiexp(chisq::Function,
                xp::Vector{Float64}, 
                data::Vector{uwreal};
                W = Vector{Float64}())

    n = length(xp)   # Number of fit parameters
    m = length(data) # Number of data

    xav = zeros(Float64, n+m)
    for i in 1:n
        xav[i] = xp[i]
    end
    for i in n+1:n+m
        xav[i] = data[i-n].mean
    end
    ccsq(x::Vector) = chisq(view(x, 1:n), view(x, n+1:n+m)) 
    cfg = ForwardDiff.HessianConfig(ccsq, xav, Chunk{8}());

    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)
        
    cse = 0.0
    if (m-n > 0)
        if (length(W) == 0)
            Ww = zeros(Float64, m)
            for i in 1:m
                if (data[i].err == 0.0)
                    uwerr(data[i])
                    if (data[i].err == 0.0)
                        error("Zero error in fit data")
                    end
                end
                Ww[i] = 1.0 / data[i].err^2
            end
        else
            Ww = W
        end

        cse = chiexp(hess, data, Ww)        
    end

    return cse
end

function fit_error(chisq::Function,
                   xp::Vector{Float64}, 
                   data::Vector{uwreal};
                   W = Vector{Float64}(), chi_exp = true)

    n = length(xp)   # Number of fit parameters
    m = length(data) # Number of data

    xav = Vector{Float64}(undef, n+m)
    for i in 1:n
        xav[i] = xp[i]
    end
    for i in n+1:n+m
        xav[i] = data[i-n].mean
    end

    ccsq(x::Vector) = chisq(x[1:n], x[n+1:n+m])
    cfg = ForwardDiff.HessianConfig(ccsq, xav, Chunk{8}());

    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)
    
    hinv = LinearAlgebra.pinv(hess[1:n,1:n])
    grad = - hinv[1:n,1:n] * hess[1:n,n+1:n+m]
    
    param = addobs(data, grad, xp)

    if (!chi_exp)
        return param
    end
    
    cse = 0.0
    if (m-n > 0)
        if (length(W) == 0)
            Ww = zeros(Float64, m)
            for i in 1:m
                if (data[i].err == 0.0)
                    uwerr(data[i])
                    if (data[i].err == 0.0)
                        error("Zero error in fit data")
                    end
                end
                Ww[i] = 1.0 / data[i].err^2
            end
        else
            Ww = W
        end
        
        cse = chiexp(hess, data, Ww)        
    end

    return param, cse
end

function int_error(fint::Function,
                   a,
                   b,
                   param::Vector{uwreal})

    # First basic integral
    xp = zeros(Float64, length(param))
    for i in 1:length(xp)
        xp[i] = param[i].mean
    end
    f(x) = fint(x, xp)
    (val, foo) = QuadGK.quadgk(f, a, b)

    # Functions that return the derivatives
    # with respect to the parameters as functions
    function fp(n::Int64)
        function fd(x)
            fax(pa) = fint(x,[xp[1:n-1];[pa];xp[n+1:end]])
            return ForwardDiff.derivative(fax, xp[n])
        end
        return fd
    end

    grad = zeros(Float64, length(param))
    for i in 1:length(param)
        (grad[i], foo) = QuadGK.quadgk(fp(i), a, b)
    end

    return addobs(param, grad, val)
end

function int_error(fint::Function,
                   a::uwreal,
                   b::Float64,
                   param::Vector{uwreal})

    # First basic integral
    xp = zeros(Float64, length(param))
    for i in 1:length(xp)
        xp[i] = param[i].mean
    end
    f(x) = fint(x, xp)
    (val, foo) = QuadGK.quadgk(f, a.mean, b)

    # Functions that return the derivatives
    # with respect to the parameters as functions
    function fp(n::Int64)
        function fd(x)
            fax(pa) = fint(x,[xp[1:n-1];[pa];xp[n+1:end]])
            return ForwardDiff.derivative(fax, xp[n])
        end
        return fd
    end

    grad = zeros(Float64, length(param)+1)
    for i in 1:length(param)
        (grad[i], foo) = QuadGK.quadgk(fp(i), a.mean, b)
    end
    grad[end] = fint(a.mean, xp)

    return addobs([param; a], grad, val)
end

function int_error(fint::Function,
                   a::uwreal,
                   b::uwreal,
                   param::Vector{uwreal})

    # First basic integral
    xp = zeros(Float64, length(param))
    for i in 1:length(xp)
        xp[i] = param[i].mean
    end
    f(x) = fint(x, xp)
    (val, foo) = QuadGK.quadgk(f, a.mean, b.mean)

    # Functions that return the derivatives
    # with respect to the parameters as functions
    function fp(n::Int64)
        function fd(x)
            fax(pa) = fint(x,[xp[1:n-1];[pa];xp[n+1:end]])
            return ForwardDiff.derivative(fax, xp[n])
        end
        return fd
    end

    grad = zeros(Float64, length(param)+2)
    for i in 1:length(param)
        (grad[i], foo) = QuadGK.quadgk(fp(i), a.mean, b.mean)
    end
    grad[end-1] = fint(a.mean, xp)
    grad[end]   = fint(b.mean, xp)

    return addobs([param; a; b], grad, val)
end

function int_error(fint::Function,
                   a::Float64,
                   b::uwreal,
                   param::Vector{uwreal})

    # First basic integral
    xp = zeros(Float64, length(param))
    for i in 1:length(xp)
        xp[i] = param[i].mean
    end
    f(x) = fint(x, xp)
    (val, foo) = QuadGK.quadgk(f, a, b.mean)

    # Functions that return the derivatives
    # with respect to the parameters as functions
    function fp(n::Int64)
        function fd(x)
            fax(pa) = fint(x,[xp[1:n-1];[pa];xp[n+1:end]])
            return ForwardDiff.derivative(fax, xp[n])
        end
        return fd
    end

    grad = zeros(Float64, length(param)+1)
    for i in 1:length(param)
        (grad[i], foo) = QuadGK.quadgk(fp(i), a, b.mean)
    end
    grad[end]   = fint(b.mean, xp)

    return addobs([param; b], grad, val)
end
