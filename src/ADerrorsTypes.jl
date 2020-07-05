###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    ADerrorsTypes.jl
### created: Wed Jun 17 13:00:32 2020
###                               


mutable struct cfdata
    var::Float64
    taui::Float64
    dtaui::Float64

    iw::Int64
    
    gamm::Vector{Float64}
    drho::Vector{Float64}
end

mutable struct uwreal
    mean::Float64
    err::Float64
    derr::Float64

    prop::Array{Bool, 1}
    der::Array{Float64, 1}

    ids::Array{Int64, 1}

    cfd::Vector{cfdata}

    uwreal(a::Float64, b::Float64, c::Float64,
           d::Array{Bool, 1}, e::Array{Float64, 1},
           f::Array{Int64, 1}, g::Vector{cfdata}) = new(a, b, c, d, e, f, g) 
    function uwreal(n::Int64)
        x = new()
        x.mean = 0.0
        x.err  = 0.0
        x.derr = 0.0
        x.prop = Vector{Bool}(undef, n)
        x.der  = Vector{Float64}(undef, n)

        x.prop .= false
        x.der  .= 0.0
        return x
    end
end

mutable struct fbd
    nd::Int64
    ibn::Int64
    
    delta::Array{Float64, 1}
    ivrep::Array{Int64, 1}
    fourier::Dict{Int64,Vector{Complex{Float64}}}
end

mutable struct wspace
    fluc::Array{fbd, 1}
    nob::Int64

    map_nob::Array{Int64, 1}    # The id that corresponds to each ob
    map_ids::Dict{Int64, Int64} # For each id, a fluc index
end


Base.convert(::Type{uwreal}, x::Float64) = uwreal(x, 0.0, 0.0, 
                                            Vector{Bool}(), Vector{Float64}(), 
                                            Vector{Int64}(), Vector{cfdata}())

uwreal(v::Float64, n::Int) = uwreal(v, 0.0, 0.0,
                                    [false for i in 1:n], [0.0 for i in 1:n],
                                    Vector{Int64}(), Vector{cfdata}())
uwreal(v::Float64, prop::Vector{Bool}, der::Vector{Float64}) = uwreal(v, 0.0, 0.0,
                                    prop, der,
                                    Vector{Int64}(), Vector{cfdata}())

function Base.show(io::IO, a::uwreal)
    
    if (length(a.prop) == 0)
        print(a.mean)
        return
    end
    
    if (length(a.cfd) > 0) 
        print(a.mean, " +/- ", a.err)
    else
        print(a.mean, " (Error not available... maybe run uwerr)")
    end
end
