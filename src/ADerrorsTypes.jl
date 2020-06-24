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
end

mutable struct fbd
    nd::Int64
    ibn::Int64
    
    delta::Array{Float64, 1}
    ivrep::Array{Int64, 1}
    is::Array{Int64, 1}
    ie::Array{Int64, 1}
    fourier::Dict{Int64,Vector{Complex{Float64}}}
end

mutable struct wspace
    fluc::Array{fbd, 1}
    nob::Int64

    map_nob::Array{Int64, 1}    # The id that corresponds to each ob
    map_ids::Dict{Int64, Int64} # For each id, a fluc index
end


uwreal(v::Float64, n::Int) = uwreal(v, 0.0, 0.0,
                                    [false for i in 1:n], [0.0 for i in 1:n],
                                    Vector{Int64}(), Vector{cfdata}())
uwreal(v::Float64, prop::Vector{Bool}, der::Vector{Float64}) = uwreal(v, 0.0, 0.0,
                                    prop, der,
                                    Vector{Int64}(), Vector{cfdata}())
