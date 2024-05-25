###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    ADerrorsTools.jl
### created: Fri Jun 26 21:48:35 2020
###                               

"""
    derivative(a::uwreal, b::uwreal)

Determines the derivate of the observable `a` with respect to observable `b` (i.e. ``\frac{{\rm d} a}{{\rm d} b}). Observable `b` must depend on only a single ensemble (typically `b` is just a primary observable). 
```@example
using ADerrors # hide

# Put some data 
a = uwreal(1.2 .+ randn(2000), "Ensemble test A", ["R0", "Rep 1", "Last rep"], [1000, 550, 450])
b = uwreal(0.8 .+ randn(2034), "Ensemble test B", [1000, 584, 450])

d = sin(a+b)
println("Derivative of d w.r.t a:      ", derivative(d, a))
println("Derivative of d w.r.t log(a): ", derivative(d, log(a)))
```
"""
function derivative(a::uwreal, p::uwreal)

    if count(p.prop .== true) != 1
        error("I do not know how to compute this")
    end
    
    for i in 1:min(length(a.der),length(p.der))
        if (p.prop[i] && a.prop[i])
            return res = a.der[i]/p.der[i]
        end 
    end 
    
    return 0.0


end


"""
    cobs(avgs::Vector{Float64}, Mcov::Array{Float64, 2}, str::String)
    cobs(avgs::Vector{Float64}, Mcov::Array{Float64, 2}, ids::Vector{Int64})

Returns a vector of `uwreal` such that their mean values are `avgs[:]` and their covariance is `Mcov[:,:]`. In order to construct these observables `n=length(avgs)` ensemble ID are used. These are generated either from the string `str` or have to be specified in the vector `ids[:]`.
```@example
using ADerrors # hide
# Put some average values and covariance 
avg = [16.26, 0.12, -0.0038]
Mcov = [0.478071 -0.176116 0.0135305
        -0.176116 0.0696489 -0.00554431
        0.0135305 -0.00554431 0.000454180]

# Produce observables with ensemble ID 
# [1, 2001, 32]. Do error analysis
p = cobs(avg, Mcov, "GF beta function parameters")
uwerr.(p)

# Check central values are ok
avg2 = value.(p)
println("Better be zero: ", sum((avg.-avg2).^2))

# Check that the covariance is ok
Mcov2 = cov(p)
println("Better be zero: ", sum((Mcov.-Mcov2).^2))
```
"""
function cobs(avgs::Vector{Float64}, cov::Array{Float64, 2}, ids::Vector{Int64})


    n = length(avgs)
    p = Vector{uwreal}(undef, n)
    try
        ch = LinearAlgebra.cholesky(cov)
        
        for j in 1:n
            p[j] = uwreal([avgs[j], ch.L[j, 1]], ids[1])
            for i in 2:n
                p[j] = p[j] + uwreal([0.0, ch.L[j,i]], ids[i])
            end
        end

        return p
    catch
        usvt = LinearAlgebra.svd(cov)
        SS = similar(usvt.S)
        for i in eachindex(SS)
            if (usvt.S[i] > 1.0e-6)
                SS[i] = usvt.S[i]
            else
                SS[i] = 0.0
            end
        end
        A    = usvt.U * LinearAlgebra.Diagonal(sqrt.(usvt.S))

        for j in 1:n
            p[j] = uwreal([avgs[j], A[j, 1]], ids[1])
            for i in 2:n
                p[j] = p[j] + uwreal([0.0, A[j,i]], ids[i])
            end
        end
    end
    
    return p
end

function cobs(avgs::Vector{Float64}, cov::Array{Float64, 2}, sids::Vector{String})


    n = length(avgs)
    p = Vector{uwreal}(undef, n)
    try
        ch = LinearAlgebra.cholesky(cov)
        
        for j in 1:n
            p[j] = uwreal([avgs[j], ch.L[j, 1]], sids[1])
            for i in 2:n
                p[j] = p[j] + uwreal([0.0, ch.L[j,i]], sids[i])
            end
        end

        return p
    catch
        usvt = LinearAlgebra.svd(cov)
        SS = similar(usvt.S)
        for i in eachindex(SS)
            if (usvt.S[i] > 1.0e-6)
                SS[i] = usvt.S[i]
            else
                SS[i] = 0.0
            end
        end
        A    = usvt.U * LinearAlgebra.Diagonal(sqrt.(usvt.S))

        for j in 1:n
            p[j] = uwreal([avgs[j], A[j, 1]], sids[1])
            for i in 2:n
                p[j] = p[j] + uwreal([0.0, A[j,i]], sids[i])
            end
        end
    end
    
    return p
end

function cobs(avgs::Vector{Float64}, cov::Array{Float64, 2}, str::String)

    n = length(avgs)
    sids = Vector{String}(undef, n)
    for i in 1:n
        sids[i] = Printf.@sprintf "%s %8.8d" str i
    end

    return cobs(avgs, cov, sids)
end


function addobs(a::Vector{uwreal}, der::Vector{Float64}, mvl::Float64)

    n = 0
    for i in 1:length(a)
        n = max(n,length(a[i].der))
    end
    p = [false for i in 1:n]
    d = zeros(Float64, n)

    for k in 1:n
        for j in 1:length(a)
            if (k <= length(a[j].prop))
                p[k] = p[k] || a[j].prop[k]
                if (a[j].prop[k])
                    d[k] = d[k] + der[j]*a[j].der[k]
                end
            end
        end
    end
    
    return uwreal(mvl, p, d)
end

function addobs(a::Vector{uwreal}, der::Array{Float64, 2}, mvl::Vector{Float64})

    n = 0
    for i in 1:length(a)
        n = max(n,length(a[i].der))
    end

    x = Vector{uwreal}(undef, size(der, 1))
    for i in 1:size(der, 1)
        p = [false for i in 1:n]
        d = zeros(Float64, n)
        for k in 1:n
            for j in 1:length(a)
                if (k <= length(a[j].prop))
                    p[k] = p[k] || a[j].prop[k]
                    if (a[j].prop[k])
                        d[k] = d[k] + der[i,j]*a[j].der[k]
                    end
                end
            end
        end
        x[i] = uwreal(mvl[i], p, d)
    end
    
    return x
end
