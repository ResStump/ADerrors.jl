###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    ADerrorsIO.jl
### created: Tue Jul 1 19:40:39 2020
###                               

function find_mcid(a::uwreal, mcid::Int64)

    if (length(a.cfd) == 0)
        return nothing
    else
        for i in 1:length(a.ids)
            if (a.ids[i] == mcid)
                return i
            end
        end
    end
    
    return nothing
end

"""
    err(a::uwreal)

Returns the error of the `uwreal` variable `a`. It is assumed that `uwerr` has been run on the variable so that an error is available. Otherwise an error message is printed.
```@example
using ADerrors # hide
a = uwreal([1.2, 0.2], 12)   # a = 1.2 +/- 0.2
uwerr(a)
println("a has error: ", err(a))
```
"""
function err(a::uwreal)
    if (length(a.cfd) == 0)
        error("No error available... maybe run uwerr")
    end
    return a.err
end

"""
    value(a::uwreal)

Returns the (mean) value of the `uwreal` variable `a`
```@example
using ADerrors # hide
a = uwreal([1.2, 0.2], 12)   # a = 1.2 +/- 0.2
uwerr(a)
println("a has central value: ", value(a))
```
"""
value(a::uwreal)            = a.mean

"""
    derror(a::uwreal)

Returns an estimate of teh error of the error of the `uwreal` variable `a`. It is assumed that `uwerr` has been run on the variable so that an error is available. Otherwise an error message is printed.
```@example
using ADerrors # hide
a = uwreal([1.2, 0.2], 12)   # a = 1.2 +/- 0.2
uwerr(a)
println("a has error of the error: ", derror(a))
```
"""
function derror(a::uwreal)
    if (length(a.cfd) == 0)
        error("No error available... maybe run uwerr")
    end
    return a.derr
end

"""
    taui(a::uwreal, mcid)

Returns the value of tauint for the ensemble `mcid`. It is assumed that `uwerr` has been run on the variable and that `mcid` contributes to the observable `a`. Otherwise an error message is printed. `mcid` can be either an `Int64` (the proper ensemble ID), or a `String` (the ensemble tag).
```@example
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

a = uwreal(x.^2, "Some simple ensemble")
uwerr(a)
println("Error analysis result: ", a, " (tauint = ", taui(a, "Some simple ensemble"), ")")
```
"""
function taui(a::uwreal,   mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == nothing)
        error("No error available... maybe run uwerr")
    else
        return a.cfd[idx].taui
    end
end

"""
    dtaui(a::uwreal, mcid)

Returns an estimate on the error of tauint for the ensemble `mcid`. It is assumed that `uwerr` has been run on the variable and that `mcid` contributes to the observable `a`. Otherwise an error message is printed.
```@example
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

a = uwreal(x.^2, "Some simple ensemble")
uwerr(a)
println("Error analysis result: ", a, 
        " (tauint = ", taui(a, "Some simple ensemble"), " +/- ", dtaui(a, "Some simple ensemble"), ")")
```
"""
function dtaui(a::uwreal,  mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == nothing)
        error("No error available... maybe run uwerr")
    else
        return a.cfd[idx].dtaui
    end
end

"""
    window(a::uwreal, mcid)

Returns the summation window for the ensemble `mcid`. It is assumed that `uwerr` has been run on the variable and that `mcid` contributes to the observable `a`. Otherwise an error message is printed.
```@example
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

a = uwreal(x.^2, "Some simple ensemble")
uwerr(a)
println("Error analysis result: ", a, 
        " (window = ", window(a, "Some simple ensemble"), ")")
```
"""
function window(a::uwreal, mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == nothing)
        error("No error available... maybe run uwerr")
    else
        return a.cfd[idx].iw
    end
end

"""
    rho(a::uwreal, mcid)

Returns the normalized autocorrelation function of `a` for the ensemble `mcid`. It is assumed that `uwerr` has been run on the variable and that `mcid` contributes to the observable `a`. Otherwise an error message is printed.
```@example
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

a = uwreal(x.^2, "Some simple ensemble")
uwerr(a)
v = rho(a, "Some simple ensemble")
for i in 1:length(v)
    println(i, " ", v[i])
end
```
"""
function rho(a::uwreal, mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == nothing)
        error("No error available... maybe run uwerr")
    else
        return a.cfd[idx].gamm ./ a.cfd[idx].gamm[1]
    end
end

"""
    drho(a::uwreal, mcid)

Returns an estimate of the error on the normalized autocorrelation function of `a` for the ensemble `mcid`. It is assumed that `uwerr` has been run on the variable and that `mcid` contributes to the observable `a`. Otherwise an error message is printed.
```@example
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

a = uwreal(x.^2, "Some simple ensemble")
uwerr(a)
v  =  rho(a, "Some simple ensemble")
dv = drho(a, "Some simple ensemble")
for i in 1:length(v)
    println(i, " ", v[i], " +/- ", dv[i])
end
```
"""
function drho(a::uwreal, mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == nothing)
        error("No error available... maybe run uwerr")
    else
        return a.cfd[idx].drho
    end
end

function read_bdio(fb, ws::wspace, mapids::Dict{Int64, Int64})

    dfoo = zeros(Float64, 1)
    ifoo = zeros(Int32, 1)
    BDIO.BDIO_read(fb, dfoo)
    BDIO.BDIO_read(fb, ifoo)

    nid::Int32 = ifoo[1]

    p = [false for i in 1:ws.nob+nid]
    d = zeros(Float64, ws.nob+nid)
    p[ws.nob+1:end] .= true
    d[ws.nob+1:end] .= 1.0

    nds  = zeros(Int32, nid)
    nrep = zeros(Int32, nid)
    ids  = zeros(Int32, nid)
    itmp = zeros(Int32, nid)

    BDIO.BDIO_read(fb, nds)
    BDIO.BDIO_read(fb, nrep)
    ivrep = zeros(Int32, sum(nrep))
    BDIO.BDIO_read(fb, ivrep)
    BDIO.BDIO_read(fb, ids)

    dfl = zeros(Float64, nid)
    BDIO.BDIO_read(fb, itmp)
    for i in 1:2
        BDIO.BDIO_read(fb, dfl)
    end

    is = 1
    for i in 1:nid
        ie = is + nrep[i] - 1
        if (sum(ivrep[is:ie]) != nds[i])
            throw("Replica sum does not match number of measurements")
        end

        dfl = zeros(Float64, nds[i])
        BDIO.BDIO_read(fb, dfl)
        id = convert(Int64, ids[i])
        add_DB(dfl, get(mapids, id, id), convert(Vector{Int64}, ivrep[is:ie]), ws, false)
        
        is = ie + 1
    end

    if BDIO.BDIO_eor(fb)
        for i in 1:nid
            id = convert(Int64, ids[i])
            add_maps(id, ws)
            get_name_from_id(id, ws)
        end
    else
        name = BDIO.BDIO_read_str(fb)

        for i in 1:nid
            BDIO.BDIO_read(fb, ifoo)
            str = BDIO.BDIO_read_str(fb)

            id = get_id_from_name(str, ws)
            add_maps(id, ws)
        end
    end
    
    return uwreal(dfoo[1], p, d)
end

function write_bdio(p::uwreal, fb, iu::Int, ws::wspace; name="NO NAME")

    BDIO.BDIO_start_record!(fb, BDIO.BDIO_BIN_GENERIC, convert(Int32, iu), true)
    nid = convert(Int32, unique_ids!(p, ws))

    BDIO.BDIO_write!(fb, [p.mean], true)
    BDIO.BDIO_write!(fb, [nid], true)
    ntv = Vector{Int32}(undef, nid)
    for j in 1:nid
        ntv[j] = convert(Int32, ws.fluc[ws.map_ids[p.ids[j]]].nd)
    end
    BDIO.BDIO_write!(fb, ntv, true)

    for j in 1:nid
        ntv[j] = convert(Int32, length(ws.fluc[ws.map_ids[p.ids[j]]].ivrep))
    end
    BDIO.BDIO_write!(fb, ntv, true)

    for j in 1:nid
        BDIO.BDIO_write!(fb, convert(Vector{Int32}, ws.fluc[ws.map_ids[p.ids[j]]].ivrep), true)
    end

    BDIO.BDIO_write!(fb, convert(Vector{Int32}, p.ids), true)

    for j in 1:nid
        ntv[j] = convert(Int32, div(maximum(ws.fluc[ws.map_ids[p.ids[j]]].ivrep), 2))
    end
    BDIO.BDIO_write!(fb, ntv, true)

    BDIO.BDIO_write!(fb, zeros(Float64, nid), true)
    BDIO.BDIO_write!(fb, [DEFAULT_STAU for n in 1:nid], true)

    for j in 1:nid
        nd = ws.fluc[ws.map_ids[p.ids[j]]].nd
        dt = zeros(Float64, nd)
        for i in 1:length(p.prop)
            if (p.prop[i] && (ws.map_nob[i] == p.ids[j]))
                dt .= dt .+ p.der[i] .* ws.fluc[i].delta
            end
        end
        BDIO.BDIO_write!(fb, dt, true)
    end

    BDIO.BDIO_write!(fb, name*"\0")
    for i in 1:nid
        BDIO.BDIO_write!(fb, [convert(Int32, p.ids[i])])
        BDIO.BDIO_write!(fb, get_name_from_id(p.ids[i], ws)*"\0")
    end
    
    BDIO.BDIO_write_hash!(fb)

    return true
end

"""
    details(a::uwreal; io::IO=stdout, names::Dict{Int64, String} = Dict{Int64, String}())

Write out a detailed information on the error of `a`.

## Arguments

Optionally one can pass as a keyword argument (`io`) the `IO` stream to write to.

## Example
```@example
using ADerrors # hide
a = uwreal(rand(2000),   "Ensemble A12")
b = uwreal([1.2, 0.023], "Ensemble XYZ")
c = uwreal([5.2, 0.03],  "Ensemble RRR")
d = a + b - c
uwerr(d)

details(d)
```
"""
function details(a::uwreal, ws::wspace, io::IO=stdout)
    
    if (length(a.prop) == 0)
        print(a.mean)
        return
    end

    if (length(a.cfd) > 0) 
        println(io, a.mean, " +/- ", a.err)
        println(io, " ## Number of error sources: ", length(a.ids))

        n = 0
        for i in 1:length(a.cfd)
            idx  = ws.map_ids[a.ids[i]]
            if (ws.fluc[idx].nd > 1)
                n = n + 1
            end
        end
        println(io, " ## Number of MC ids       : ", n)
        println(io, " ## Contribution to error  :               Ensemble  [%]     [MC length]")

        
        truncate_ascii(s::String,n::Int) = s[1:min(sizeof(s),n)]
        ntrunc = 45
        v = zeros(length(a.cfd))
        for i in eachindex(v)
            v[i] = a.cfd[i].var
        end
        ip = sortperm(v, rev=true)
        for i in 1:length(a.cfd)
            idx  = ws.map_ids[a.ids[ip[i]]]
            sndt = join(ws.fluc[idx].ivrep, ",")
            sid  = truncate_ascii(get_name_from_id(a.ids[ip[i]], ws), ntrunc)
            if (ws.fluc[idx].nd > 1)
                Printf.@printf(io, "  #  %45s %6.2f   %s\n",
                        sid, 100.0 .* a.cfd[ip[i]].var ./ a.err^2, sndt)
            else
                Printf.@printf(io, "  #  %45s %6.2f            -\n",
                        sid, 100.0 .* a.cfd[ip[i]].var ./ a.err^2)
            end
        end
    else
        print(io, a.mean, " (Error not available... maybe run uwerr)")
    end
end

details(a::uwreal; io::IO=stdout) = details(a, wsg, io)

"""
    read_uwreal(fb[, map_ids::Dict{Int64, Int64}])

Given a `BDIO` file handler `fb`, this routine returns the observable stored in the current record. Optionally, ensemble ID stored in the file can be changed at read time by passing a dictionary `map_ids`
```@example
using ADerrors # hide
using BDIO
a = uwreal(rand(2000), 12)

fb = BDIO_open("/tmp/foo.bdio", "w", "Test file")
write(a, fb, 8)
BDIO_close(fb)

# Open the file and move to first record
fb = BDIO_open("/tmp/foo.bdio", "r")
BDIO_seek!(fb)

# Read observable
changeID_12_by_120 = Dict(12 => 120)
b = read_uwreal(fb, changeID_12_by_120)

# Check
c = a - b
uwerr(c)
print("Better be zero: ")
details(c)
```
"""
read_uwreal(fb, mapids::Dict{Int64, Int64} = Dict{Int64, Int64}())  = read_bdio(fb, ADerrors.wsg, mapids)

"""
    write_uwreal(p::uwreal, fb, iu::Int)

Given a `BDIO` file handler `fb`, this writes the observable `p` in a BDIO record with user info `iu`.
```@example
using ADerrors # hide
using BDIO
a = uwreal(rand(2000), 12)

# Create a BDIO file and write observable with user info 8.
fb = BDIO_open("/tmp/foo.bdio", "w", "Test file")
write(a, fb, 8)
BDIO_close(fb)

# Open the file and move to first record
fb = BDIO_open("/tmp/foo.bdio", "r")
BDIO_seek!(fb)

# Read observable
b = read_uwreal(fb)

# Check
c = a - b
uwerr(c)
println("Better be zero: ", c)
```
"""
write_uwreal(p::uwreal, fb, iu::Int) = write_bdio(p::uwreal, fb, iu::Int, ADerrors.wsg)

