
a = uwreal([1.0, 0.1], 1)
b = uwreal([0.5, 0.023], 23)

c = 1.0 + sin(a+b)
d = sin(a)*cos(b) + cos(a)*sin(b) - 3.0

let e = c-d, nmax = 1000
    for i in 1:nmax
        e = e + c-d
    end
    uwerr(e)
    println(e.mean, " +/- ", e.err)
    ( (abs(err(e)) < 1.0E-10) && (abs(value(e)-4.0*(nmax+1.0)) < 1.0E-10) )
end
