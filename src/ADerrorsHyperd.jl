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

const LG2  = 0.6931471805599453094172321214581765680755001343602552
const LG10 = 2.3025850929940456840179914546843642076011014886287729

function Base.sqrt(h::hyperd)
    v  = sqrt(h.v)
    d1 = 0.5/v
    dd = -0.25/v^3
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.log(h::hyperd)
    v  = log(h.v)
    d1 = 1.0/h.v
    dd = -d1/h.v
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.log2(h::hyperd)
    v  = log2(h.v)
    d1 = 1.0/(h.v*LG2)
    dd = -d1/h.v
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.log10(h::hyperd)
    v  = log2(h.v)
    d1 = 1.0/(h.v*LG10)
    dd = -d1/h.v
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.exp(h::hyperd)
    v  = exp(h.v)
    d1 = v
    dd = v
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.exp2(h::hyperd)
    v  = exp2(h.v)
    d1 = v * LG2
    dd = d1 * LG2
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.exp10(h::hyperd)
    v  = exp10(h.v)
    d1 = v * LG10
    dd = d1 * LG10
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.sin(h::hyperd)
    v  = sin(h.v)
    d1 = cos(h.v)
    dd = -v
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.cos(h::hyperd)
    v  =  cos(h.v)
    d1 = -sin(h.v)
    dd = -v
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.tan(h::hyperd)
    v  = tan(h.v)
    d1 = 1.0 + v*v
    dd = 2.0 * v*d1
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.sec(h::hyperd)
    v  = sec(h.v)
    d1 = v*tan(h.v)
    dd = v*(tan(h.v)^2 + v^2)
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.csc(h::hyperd)
    v  = csc(h.v)
    d1 = -v * cot(h.v)
    dd = v*(cot(h.v)^2 + v^2)
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end

function Base.cot(h::hyperd)
    v  = cot(h.v)
    d1 =  -(1 + v^2)
    dd = -2.0*v*d1
    return hyperd(v, d1*h.d1, d1*h.d2, dd*h.d1*h.d2 + d1*h.dd)
end



