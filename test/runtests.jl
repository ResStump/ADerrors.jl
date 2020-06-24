#/usr/bin/env julia14

# Start test script
using ADerrors
using Test

println("Test [test1.jl]")
@time @test include("test1.jl")

