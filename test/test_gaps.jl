using ADerrors, Plots
pgfplotsx();

# Generate some correlated data
eta  = randn(10000);
x    = Vector{Float64}(undef, 10000);
x[1] = 0.0;
for i in 2:10000
    x[i] = x[i-1] + 0.2*eta[i]
    if abs(x[i]) > 1.0
        x[i] = x[i-1]
    end
end

# x^2 only measured on odd configurations
x2 = uwreal(x[1:2:9999].^2,  "Random walk in [-1,1]", collect(1:2:9999),  10000)
# x^4 only measured on even configurations
x4 = uwreal(x[2:2:10000].^4, "Random walk in [-1,1]", collect(2:2:10000), 10000)

rat = x2/x4 # This is perfectly fine
try
    uwerr(rat); # Use automatic window: might fail!
catch e
    println("Automatic window fails")
end
    
iw = window(rat, "Random walk in [-1,1]")
r  = rho(rat, "Random walk in [-1,1]");
dr = drho(rat, "Random walk in [-1,1]");
plot(collect(1:100), 
	r[1:100], 
	yerr = dr[1:100], 
	seriestype = :scatter, title = "Chosen Window: " * string(iw), label="autoCF")
savefig("rat_cf.png") # hide

wpm = Dict{String, Vector{Float64}}()
wpm["Random walk in [-1,1]"] = [50.0, -1.0, -1.0, -1.0]
uwerr(rat, wpm) # repeat error analysis with our choice (window=50)
println("Ratio:   ", rat)

prod = x2*x4
uwerr(prod)
println("Product: ", prod)
iw = window(prod, "Random walk in [-1,1]")
r  = rho(prod, "Random walk in [-1,1]");
dr = drho(prod, "Random walk in [-1,1]");
plot(collect(1:2*iw), 
	r[1:2*iw], 
	yerr = dr[1:2*iw], 
	seriestype = :scatter, title = "Chosen Window: " * string(iw), label="autoCF")
savefig("prod_cf.png") # hide
