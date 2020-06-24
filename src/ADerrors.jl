###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    ADerrors.jl
### created: Wed Jun 17 13:00:26 2020
###                               

module ADerrors

import ForwardDiff, Statistics, FFTW#, BDIO

# Include data types
include("ADerrorsTypes.jl")

# Include computation of autoCF
include("ADerrorsCF.jl")

# Math operations
include("ADerrorsMath.jl")

# I/O
include("ADerrorsIO.jl")
export err, value, derror, taui, dtaui, window
export uwreal, uwerr

end # module
