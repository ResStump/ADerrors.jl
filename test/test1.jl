using ADerrors # hide
# Generate some correlated data
eta  = randn(1000)
x    = Vector{Float64}(undef, 1000)
x[1] = 0.0
for i in 2:1000
    x[i] = x[i-1] + eta[i]
    if abs(x[i]) > 1.0
        x[i] = x[i-1]
    end
end

# Load the data in a uwreal
a = uwreal(x.^2, "Random walk in [-1,1]")
wpm = Dict{String,Vector{Float64}}()

# Use default analysis (stau = 4.0)
uwerr(a)
println("default:                   ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")

# This will still do default analysis because 
# a does not depend on emsemble foo
wpm["Ensemble foo"] = [-1.0, 8.0, -1.0, 145.0]
uwerr(a, wpm)
println("default:                   ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")

# Fix the summation window to 1 (i.e. uncorrelated data)
wpm["Random walk in [-1,1]"] = [1.0, -1.0, -1.0, -1.0]
uwerr(a, wpm)
println("uncorrelated:              ",  a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")

# Use stau = 1.5
wpm["Random walk in [-1,1]"] = [-1.0, 1.5, -1.0, -1.0]
uwerr(a, wpm)
println("stau = 1.5:                ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")

# Use fixed window 15 and add tail with texp = 100.0
wpm["Random walk in [-1,1]"] = [15.0, -1.0, -1.0, 100.0]
uwerr(a, wpm)
println("Fixed window 15, texp=100: ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")

# Sum up to the point that the signal in Gamma is 
# 1.5 times the error and add a tail with texp = 10.0
wpm["Random walk in [-1,1]"] = [-1.0, -1.0, 1.5, 30.0]
uwerr(a, wpm)
println("signal/noise=1.5, texp=10: ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")

(0 == 0)
