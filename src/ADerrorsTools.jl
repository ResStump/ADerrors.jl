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
    cobs(avgs::Vector{Float64}, Mcov::Array{Float64, 2}, ids::Vector{Int64})

Returns a vector of `uwreal` such that their mean values are `avgs[:]` and their covariance is `Mcov[:,:]`. In order to construct these observables `n=length(avgs)` ensemble ID are used. These have to be specified in the vector `ids[:]`.
```@example
using ADerrors # hide
# Put some average values and covariance 
avg = [16.26, 0.12, -0.0038]
Mcov = [0.478071 -0.176116 0.0135305
        -0.176116 0.0696489 -0.00554431
        0.0135305 -0.00554431 0.000454180]

# Produce observables with ensemble ID 
# [1, 2001, 32]. Do error analysis
p = cobs(avg, Mcov, [1, 2001, 32])
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

    ch = LinearAlgebra.cholesky(cov)

    n = length(avgs)
    p = Vector{uwreal}(undef, n)
    for j in 1:n
        p[j] = uwreal([avgs[j], ch.L[j, 1]], ids[1])
        for i in 2:n
            p[j] = p[j] + uwreal([0.0, ch.L[j,i]], ids[i])
        end
    end
    return p
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
