###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    ADerrorsMath.jl
### created: Thu Jun 18 08:13:30 2020
###                               

for op in (:sin, :cos, :tan, :log, :exp, :sqrt, :sind, :cosd, :tand, :sinpi, :cospi, :sinh, :cosh, :tanh, :asin, :acos, :atan, :asind, :acosd, :atand, :sec, :csc, :cot, :secd, :cscd, :cotd, :asec, :acsc, :acot, :asecd, :acscd, :acotd, :sech, :csch, :coth, :asinh, :acosh, :atanh, :asech, :acsch, :acoth, :sinc, :cosc, :deg2rad, :rad2deg, :log2, :log10, :log1p, :exp2, :exp10, :expm1, :-)
    @eval function Base.$op(a::uwreal)
        
        return uwreal(Base.$op(a.mean),
                      a.prop, ForwardDiff.derivative($op, a.mean)*a.der)
    end
end 
Base.:+(a::uwreal) = a

for op in (:+, :-, :*, :/, :^, :atan, :hypot)
    @eval function Base.$op(a::uwreal, b::uwreal)

        function fvec(x::Vector)
            return Base.$op(x[1], x[2])
        end
        x = [a.mean, b.mean]
        cfg  = GradientConfig(fvec, x, Chunk{2}());
        grad = ForwardDiff.gradient(fvec, x, cfg)
        
        if (length(a.der) > length(b.der))
            p = similar(a.prop)
            d = similar(a.der)
            @inbounds for i in 1:length(b.der)
                d[i] = grad[1]*a.der[i] + grad[2]*b.der[i]
                p[i] = a.prop[i] || b.prop[i]
            end
            @inbounds for i in length(b.der)+1:length(a.der)
                d[i] = grad[1]*a.der[i]
                p[i] = a.prop[i]
            end
        else
            p = similar(b.prop)
            d = similar(b.der)
            @inbounds for i in 1:length(a.der)
                d[i] = grad[1]*a.der[i] + grad[2]*b.der[i]
                p[i] = a.prop[i] || b.prop[i]
            end
            @inbounds for i in length(a.der)+1:length(b.der)
                d[i] = grad[2]*b.der[i]
                p[i] = b.prop[i]
            end
        end

        return uwreal(Base.$op(a.mean,b.mean), p, d)
    end

    @eval function Base.$op(a::uwreal, b::Number)

        function fvec(x::Vector)
            return Base.$op(x[1], b)
        end
        grad = ForwardDiff.gradient(fvec,[a.mean])
        
        return uwreal(Base.$op(a.mean,b), a.prop, grad[1]*a.der)
    end

    @eval function Base.$op(b::Number, a::uwreal)

        function fvec(x::Vector)
            return Base.$op(b, x[1])
        end
        grad = ForwardDiff.gradient(fvec,[a.mean])
        
        return uwreal(Base.$op(b,a.mean), a.prop, grad[1]*a.der)
    end
end

for op in (:<, :>, :≤, :≥, :≠)
    @eval Base.$op(a::uwreal, b::uwreal) = Base.$op(a.mean, b.mean)
    @eval Base.$op(a::uwreal, b::Number) = Base.$op(a.mean, b)
    @eval Base.$op(a::Number, b::uwreal) = Base.$op(a, b.mean)
end

