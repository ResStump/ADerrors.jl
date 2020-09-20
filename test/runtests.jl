#/usr/bin/env julia14

# Start test script
using ADerrors
using Test

println("Test [test1.jl]")
@time @test include("test1.jl")

println("Test [test2.jl]")
@time @test include("test2.jl")

println("Test [test_cov1.jl]")
@time @test include("test_cov1.jl")

println("Test [test_cov2.jl]")
@time @test include("test_cov2.jl")

println("Test [test_trcov.jl]")
@time @test include("test_trcov.jl")

