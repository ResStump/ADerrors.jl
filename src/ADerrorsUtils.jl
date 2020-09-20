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
data[1] = uwreal([1.0, 0.2],   "Var A")
data[2] = uwreal([1.2, 0.023], "Var B")
data[3] = uwreal(rand(1000),   "White noise ensemble")

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
print("Better be zero (with zero error): ")
details(z)
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

function chiexp(hess::Array{Float64, 2}, data::Vector{uwreal}, W::Vector{Float64}, wpm::Dict{Int64,Vector{Float64}})

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

    return trcov(Px, data, wpm)
end
    
function chiexp(hess::Array{Float64, 2}, data::Vector{uwreal}, W::Array{Float64, 2}, wpm::Dict{Int64,Vector{Float64}})

    m = length(data)
    n = size(hess, 1) - m
    
    Lm = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(W))
    Li = LinearAlgebra.inv(Lm.L)
    
    hm = view(hess, 1:n, n+1:n+m)
    sm = hm * Li'
    
    maux = sm * sm'
    hi   = LinearAlgebra.pinv(maux)
    Px   = W - hm' * hi * hm

    return trcov(Px, data, wpm)
end

@doc raw"""
    chiexp(chisq::Function,
                xp::Vector{Float64}, 
                data::Vector{uwreal};
                W = Vector{Float64}())

Given a ``\chi^2(p, d)``, function of the fit parameters `p[:]` and the data `d[:]`, compute the expected value of the ``\chi^2(p, d)``.

#### Arguments

- `chisq`: Must be a function of two vectors (i.e. `chisq(p::Vector, d::Vector)`). The function is assumed to have the form

``\chi^2(p, d) = \sum_i [d_i - f_i(p)]W_{ij}[d_j - f_j(p)]``

where the function ``f_i(p)`` is an arbitrary function of the fit parameters. In simple words, the function is assumed to be quadratic in the data.
- `xp`: A vector of `Float64`. The value of the fit parameters at the minima.
- `data`: A vector of `uwreal`. The data whose fluctuations enter in the evaluation of the `chisq`.
- `W`: A matrix. The weights that enter in the evaluation of the `chisq` function. If a vector is passed, the matrix is assumed to be diagonal (i.e. **uncorrelated** fit). If no weights are passed, the routines assumes that `W` is diagonal with entries given by the inverse errors squared of the data (i.w. the `chisq` is weighted with the errors of the data). 

#### Example 
```@example
using ADerrors, Distributions # hide

# Generate correlated samples with average 0.1
npt = 12
sig = zeros(npt, npt)
dx  = zeros(npt)
for i in 1:npt
    dx[i]    = 0.01*i
    sig[i,i] = dx[i]^2
    for j in i+1:npt
        sig[i,j] = 0.0001 - 0.000005*abs(i-j)
        sig[j,i] = 0.0001 - 0.000005*abs(i-j)
    end
end
dmv = MvNormal([0.1 for n in 1:npt], sig)
vs  = rand(dmv, 1)

# Create the uwreal data that we want to 
# fit to a constant
dt = cobs(vs[:,1], sig, [100+n for n in 1:npt])

# Define the chi^2
chisq(p, d) = sum( (d .- p[1]) .^ 2 ./ dx .^2 )

# The result of an uncorrelated fit to a 
# constant is the weighted average
xp = [sum(value.(dt) ./ dx)/sum(1.0 ./ dx)]

# Compare chi^2 and expected chi^2
println("chi^2 / chi_exp^2: ", chisq(xp, value.(dt)), " / ", chiexp(chisq, xp, dt))
```
"""
function chiexp(chisq::Function,
                xp::Vector{Float64}, 
                data::Vector{uwreal},
                wpm::Dict{Int64,Vector{Float64}};
                W::Vector{Float64} = Vector{Float64}())

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
    if (n+m < 4)
        cfg = ForwardDiff.HessianConfig(ccsq, xav, Chunk{1}());
    else
        cfg = ForwardDiff.HessianConfig(ccsq, xav, Chunk{4}());
    end
        
    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)
        
    cse = 0.0
    if (m-n > 0)
        if (length(W) == 0)
            Ww = zeros(Float64, m)
            for i in 1:m
                if (data[i].err == 0.0)
                    uwerr(data[i], wpm)
                    if (data[i].err == 0.0)
                        error("Zero error in fit data")
                    end
                end
                Ww[i] = 1.0 / data[i].err^2
            end
        else
            Ww = W
        end

        cse = chiexp(hess, data, Ww, wpm)
    end

    return cse
end
chiexp(chisq::Function,
       xp::Vector{Float64}, 
       data::Vector{uwreal};
       W::Vector{Float64} = Vector{Float64}()) = 
           chiexp(chisq, xp, data, Dict{Int64,Vector{Float64}}(), W=W)
chiexp(chisq::Function,
       xp::Vector{Float64}, 
       data::Vector{uwreal},
       wpm::Dict{String,Vector{Float64}};
       W::Vector{Float64} = Vector{Float64}()) = 
           chiexp(chisq, xp, data, dict_names_to_id(wpm), W=W)



@doc raw"""

    fit_error(chisq::Function, xp::Vector{Float64}, data::Vector{uwreal}[, wpm];
                   W = Vector{Float64}(), chi_exp = true)

Given a ``\chi^2(p, d)``, function of the fit parameters `p[:]` and the data `d[:]`, this routine return the fit parameters as `uwreal` type and optionally, the expected value of ``\chi^2(p, d)``.

#### Arguments

- `chisq`: Must be a function of two vectors (i.e. `chisq(p::Vector, d::Vector)`). To determine the fit parameters, the function can be arbitrary, but for the determination of the expected ``\chi^2(p, d)`` is assumed to have the form

``\chi^2(p, d) = \sum_{ij} [d_i - f_i(p)]W_{ij}[d_j - f_j(p)]``

where the function ``f_i(p)`` is an arbitrary function of the fit parameters. In simple words, the expected ``\chi^2(p, d)`` is determined assuming that the function ``\chi^2(p, d)`` is quadratic in the data.
- `xp`: A vector of `Float64`. The value of the fit parameters at the minima.
- `data`: A vector of `uwreal`. The data whose fluctuations enter in the evaluation of the `chisq`.
- `wpm`: `Dict{Int64,Vector{Float64}}` or `Dict{String,Vector{Float64}}`. The criteria to determine the summation window. See the documentation on `uwerr` function for more details.
- `W`: A matrix. The weights that enter in the evaluation of the `chisq` function. If a vector is passed, the matrix is assumed to be diagonal (i.e. **uncorrelated** fit). If no weights are passed, the routines assumes that `W` is diagonal with entries given by the inverse errors squared of the data (i.e. the `chisq` is weighted with the errors of the data). 
- `chi_exp`: Bool type. If false, do not compute the expected ``\chi^2(p, d)``.

#### Example 
```@example
using ADerrors, Distributions # hide

# Generate correlated samples with average 0.1
npt = 12
sig = zeros(npt, npt)
dx  = zeros(npt)
for i in 1:npt
    dx[i]    = 0.01*i
    sig[i,i] = dx[i]^2
    for j in i+1:npt
        sig[i,j] = 0.0001 - 0.000005*abs(i-j)
        sig[j,i] = 0.0001 - 0.000005*abs(i-j)
    end
end
dmv = MvNormal([0.1 for n in 1:npt], sig)
vs  = rand(dmv, 1)

# Create the uwreal data that we want to 
# fit to a constant
dt = cobs(vs[:,1], sig, "Data points")

# Define the chi^2
chisq(p, d) = sum( (d .- p[1]) .^ 2 ./ dx .^2 )

# The result of an uncorrelated fit to a 
# constant is the weighted average
xp = [sum(value.(dt) ./ dx)/sum(1.0 ./ dx)]


# Propagate errors to the fit parameters and
# determine the expected chi^2
(fitp, csqexp) = fit_error(chisq, xp, dt)

uwerr.(fitp)
println(" *** FIT RESULTS ***")
print("Fit parameter:     ")
details.(fitp)
println("chi^2 / chi_exp^2: ", chisq(xp, value.(dt)), " / ", csqexp, "  (dof: ", npt-1, ")")
```
"""
function fit_error(chisq::Function,
                   xp::Vector{Float64}, 
                   data::Vector{uwreal},
                   wpm::Dict{Int64,Vector{Float64}};
                   W::Vector{Float64} = Vector{Float64}(),
                   chi_exp::Bool = true)

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
    if (n+m < 4)
        cfg = ForwardDiff.HessianConfig(ccsq, xav, Chunk{1}());
    else
        cfg = ForwardDiff.HessianConfig(ccsq, xav, Chunk{4}());
    end

    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)
    
    hm = view(hess, 1:n, 1:n)
    sm = view(hess, 1:n, n+1:n+m)
    hinv = LinearAlgebra.pinv(hm)
    grad = - hinv * sm
    
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
                    uwerr(data[i], wpm)
                    if (data[i].err == 0.0)
                        error("Zero error in fit data")
                    end
                end
                Ww[i] = 1.0 / data[i].err^2
            end
        else
            Ww = W
        end
        
        cse = chiexp(hess, data, Ww, wpm)
    end

    return param, cse
end

function fit_error(chisq::Function,
                   xp::Vector{Float64}, 
                   data::Vector{uwreal},
                   wpm::Dict{Int64,Vector{Float64}},
                   W::Array{Float64, 2};
                   chi_exp::Bool = true)

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
    if (n+m < 4)
        cfg = ForwardDiff.HessianConfig(ccsq, xav, Chunk{1}());
    else
        cfg = ForwardDiff.HessianConfig(ccsq, xav, Chunk{4}());
    end

    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)
    
    hm = view(hess, 1:n, 1:n)
    sm = view(hess, 1:n, n+1:n+m)
    hinv = LinearAlgebra.pinv(hm)
    grad = - hinv * sm
    
    param = addobs(data, grad, xp)

    if (!chi_exp)
        return param
    end
    
    cse = 0.0
    if (m-n > 0)
        cse = chiexp(hess, data, W, wpm)
    end

    return param, cse
end

fit_error(chisq::Function,
          xp::Vector{Float64}, 
          data::Vector{uwreal};
          W::Vector{Float64} = Vector{Float64}(),
          chi_exp::Bool = true) = 
              fit_error(chisq, xp, data, Dict{Int64,Vector{Float64}}(), W=W, chi_exp=chi_exp)
fit_error(chisq::Function,
          xp::Vector{Float64}, 
          data::Vector{uwreal},
          W::Array{Float64,2};
          chi_exp::Bool = true) = 
              fit_error(chisq, xp, data, Dict{Int64,Vector{Float64}}(), W=W, chi_exp=chi_exp)
fit_error(chisq::Function,
          xp::Vector{Float64}, 
          data::Vector{uwreal},
          wpm::Dict{String,Vector{Float64}};
          W::Vector{Float64} = Vector{Float64}(),
          chi_exp::Bool = true) = 
              fit_error(chisq, xp, data, dict_names_to_id(wpm), W=W, chi_exp=chi_exp)
fit_error(chisq::Function,
          xp::Vector{Float64}, 
          data::Vector{uwreal},
          wpm::Dict{String,Vector{Float64}},
          W::Array{Float64,2};
          chi_exp::Bool = true) = 
              fit_error(chisq, xp, data, dict_names_to_id(wpm), W=W, chi_exp=chi_exp)


@doc raw"""
    int_error(fint::Function, a, b, p::Vector{uwreal})

Computes the integral

``\int_a^b {\rm fint}(x; p)\, {\rm d}x``

The errors in the parameters `p[:]` and (optionally) in the limits of the intgeral (`a`, `b`) are propagated to an error of the integral. The result is returned as an `uwreal`.
```@example
using ADerrors

# Average values and covariance from
# https://inspirehep.net/literature/1477411
# here we try to reproduce the result
# Scale ratio = 21.86(42) Eq.5.6
avg = [16.26, 0.12, -0.0038]
Mcov = [ 0.478071  -0.176116   0.0135305
        -0.176116   0.0696489 -0.00554431
         0.0135305 -0.00554431 0.000454180]

p = cobs(avg, Mcov, "Beta function fit parameters")
g1s = uwreal([2.6723, 0.0064], 4)
g2s = 11.31

fint(x, p) = - (p[1] + p[2]*x^2 + p[3]*x^4)/x^3
g1 = sqrt(g1s)
g2 = sqrt(g2s)
sint = 2.0*exp(-int_error(fint, g1, g2, p))
uwerr(sint)
print("  From integral evaluation: ")
details(sint)
```
"""
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
