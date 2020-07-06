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
