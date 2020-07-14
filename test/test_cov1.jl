using ADerrors, QuadGK, ForwardDiff

# Average values and covariance from
# https://inspirehep.net/literature/1477411
# here we try to reproduce the result
# Scale ratio = 21.86(42) Eq.5.6
avg = [16.26, 0.12, -0.0038]
Mcov = [0.478071 -0.176116 0.0135305
        -0.176116 0.0696489 -0.00554431
        0.0135305 -0.00554431 0.000454180]

p = cobs(avg, Mcov, [1, 2001, 32])
g1s = uwreal([2.6723, 0.0064], 4)
g2s = 11.31

fs(a, b, p) = -p[1]/2.0 * (1.0/b-1.0/a) +
    p[2]/2.0 * log(b/a) +
    p[3]/2.0 * (b - a) 
srat = 2.0*exp(fs(g1s, g2s, p))
uwerr(srat)
println("Computation of scale factor (should be 21.86(42))")
println("  From direct evaluation:   ", srat)

fint(x, p) = - (p[1] + p[2]*x^2 + p[3]*x^4)/x^3
g1 = sqrt(g1s)
g2 = sqrt(g2s)
sint = 2.0*exp(-int_error(fint, g1, g2, p))
uwerr(sint)
print("  From integral evaluation: ")
details(sint)
( (abs(err(srat)-err(sint)) < 1.0E-10) && (abs(value(srat)-value(sint)) < 1.0E-10) )


v = value.(p)
ff(x) = fint(x, v)

a = value(g1)
b = g2
(di,foo) = quadgk(ff, a, b)
println(di)

function ftest(p)
    ff(x) = fint(x, p)
    (di,foo) = quadgk(ff, a, b)
    return di
end

function simps(f::Function, a, b, tol::Float64 = 1.0E-8)
    function trap(al, bl, f::Function, sum, n)
        nt = 2^(n-2)
        hh = (bl-al)/nt
        
        x = al+hh/2.0
        val = al + f(x)
        for i in 2:nt
            x    = x + hh
            val += f(x)
        end
        return ( sum + (bl-al)*val/nt ) / 2.0
    end
    al = min(a,b)
    bl = max(a,b)
    s1 = (bl-al)*(f(al) + f(bl))/2.0
    s  = s1
    for i in 2:1000
        s2 = trap(al, bl, f, s1, i)
        err = s - (4.0*s2 - s1)/3.0
        s   = s - err
        println(s1, " ", s2, " ", err)
        if (abs(err)<tol)
            if (b<a)
                return -s
            else
                return s
            end
        end
        s1 = s2
    end
end

function ftest2(p)
    ff(x) = fint(x, p)
    di = simps(ff, p[4], p[5])
    return di
end

vv = copy(v)
push!(vv, a)
push!(vv, b)
println(ftest(vv))
println(ftest2(vv))
hess = ForwardDiff.hessian(ftest, v)
println(hess)
hess = ForwardDiff.hessian(ftest2, vv)
println(hess)
