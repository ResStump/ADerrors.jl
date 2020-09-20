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

(abs(value(z)) < 1.0E-10 && abs(err(z)) < 1.0E-10)
