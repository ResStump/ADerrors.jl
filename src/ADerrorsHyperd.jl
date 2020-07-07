###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    ADerrorsHyperd.jl
### created: Mon Jul  6 19:58:53 2020
###                               

for op in (:sin, :cos, :tan, :log, :exp, :sqrt, :sind, :cosd, :tand, :sinpi, :cospi, :sinh, :cosh, :tanh, :asin, :acos, :atan, :asind, :acosd, :atand, :sec, :csc, :cot, :secd, :cscd, :cotd, :asec, :acsc, :acot, :asecd, :acscd, :acotd, :sech, :csch, :coth, :asinh, :acosh, :atanh, :asech, :acsch, :acoth, :sinc, :cosc, :deg2rad, :rad2deg, :log2, :log10, :log1p, :exp2, :exp10, :expm1, :-)
    @eval function Base.$op(h::hyperd)
        fvec(x::Vector) = Base.$op(x[1])
        v  = Base.$op(h.v)
        d1 = ForwardDiff.derivative($op, h.v)
        v2 = ForwardDiff.hessian(fvec, [h.v])
        return hyperd(v, d1*h.d1, d1*h.d2, v2[1]*h.d1*h.d2 + d1*h.dd)
    end
end 

Base.:+(h1::hyperd, h2::hyperd) = hyperd(h1.v+h2.v, h1.d1+h2.d1, h1.d2+h2.d2, h1.dd+h2.dd)
Base.:+(h1::hyperd, h2::Number) = hyperd(h1.v+h2, h1.d1, h1.d2, h1.dd)
Base.:+(h1::Number, h2::hyperd) = hyperd(h1+h2.v, h2.d1, h2.d2, h2.dd)

Base.:-(h1::hyperd, h2::hyperd) = hyperd(h1.v-h2.v, h1.d1-h2.d1, h1.d2-h2.d2, h1.dd-h2.dd)
Base.:-(h1::hyperd, h2::Number) = hyperd(h1.v-h2, h1.d1, h1.d2, h1.dd)
Base.:-(h1::Number, h2::hyperd) = hyperd(h1-h2.v, h2.d1, h2.d2, h2.dd)

Base.:*(h1::hyperd, h2::hyperd) = hyperd(h1.v*h2.v,
                                   h1.v*h2.d1+h1.d1*h2.v,
                                   h1.v*h2.d2+h1.d2*h2.v,
                                   h1.v*h2.dd+h1.dd*h2.v+h1.d2*h2.d1+h1.d1*h2.d2)
Base.:*(h1::hyperd, h2::Number) = hyperd(h1.v*h2, h2*h1.d1, h2*h1.d2, h2*h1.dd)
Base.:*(h1::Number, h2::hyperd) = hyperd(h1*h2.v, h1*h2.d1, h1*h2.d2, h1*h2.dd)

Base.:/(h1::hyperd, h2::hyperd) = hyperd(h1.v/h2.v,
                                   h1.d1/h2.v - h2.d1*h1.v/h2.v^2,
                                   h1.d2/h2.v - h2.d2*h1.v/h2.v^2,
                                   h1.dd/h2.v - h2.d1*h1.d2/h2.v^2 - 
                                   h2.d2*h1.d1/h2.v^2 + 
                                   h1.v*(2.0*h2.d1*h2.d2/h2.v - h2.dd)/h2.v^2)
Base.:/(h1::hyperd, h2::Number) = hyperd(h1.v/h2, h1.d1/h2, h1.d2/h2, h1.dd/h2)
Base.:/(h1::Number, h2::hyperd) = hyperd(h1/h2.v, - h2.d1*h1.v/h2.v^2, - h2.d2*h1.v/h2.v^2,
                                   h1*(2.0*h2.d1*h2.d2/h2.v - h2.dd)/h2.v^2)

Base.:^(h1::hyperd, h2::hyperd) = exp(h2*log(h1))
Base.:^(h1::hyperd, h2::Number) = exp(h2*log(h1))
function Base.:^(h1::hyperd, n::Integer)
    
    v  = h1.v^n
    if (n == 0)
        d1 = 0.0
        dd = 0.0
    elseif (n == 1)
        d1 = 1.0
        dd = 0.0
    else
        d1 = n*h1.v^(n-1)
        dd = n*(n-1)*h1.v^(n-2)
    end

    return hyperd(v, d1*h1.d1, d1*h1.d2, dd*h1.d1*h1.d2 + d1*h1.dd)
end
    
Base.:^(h1::Number, h2::hyperd) = exp(h2*log(h1))

# Missing atan, hypot


Base.zero(::Type{hyperd})          = hyperd(0.0, 0.0, 0.0, 0.0)
Base.one(::Type{hyperd})           = hyperd(1.0, 0.0, 0.0, 0.0)
Base.length(x::hyperd)             = 1
Base.iterate(x::hyperd)            = (x, nothing)
Base.iterate(x::hyperd, ::Nothing) = nothing

function hyperd_hessian!(hess::Array{Float64, 2}, f::Function, x::Vector{Float64})

    n = length(x)
    h = Vector{hyperd}(undef, n)
    for i in 1:n
        h[i] = hyperd(x[i], 0.0, 0.0, 0.0)
    end

    for i in 1:n
        h[i].d1 = 1.0
        for j in i:n
            h[j].d2 = 1.0
            res = f(h)
            hess[i,j] = res.dd
            if (j > i)
                hess[j,i] = hess[i,j]
            end
            h[j].d2 = 0.0
        end
        h[i].d1 = 0.0
    end
    return nothing
end
