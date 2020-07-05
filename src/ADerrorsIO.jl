
function find_mcid(a::uwreal, mcid::Int64)

    if (length(a.cfd) == 0)
        return 0
    else
        for i in 1:length(a.ids)
            if (a.ids[i] == mcid)
                return i
            end
        end
    end
    
    return 0
end

err(a::uwreal)              = a.err
value(a::uwreal)            = a.mean
derror(a::uwreal)           = a.derr

function taui(a::uwreal,   mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == 0)
        return 0.5
    else
        return a.cfd[idx].taui
    end
end

function dtaui(a::uwreal,  mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == 0)
        return 0.0
    else
        return a.cfd[idx].dtaui
    end
end

function window(a::uwreal, mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == 0)
        return 0
    else
        return a.cfd[idx].iw
    end
end

function rho(a::uwreal, mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == 0)
        return [1.0]
    else
        return a.cfd[idx].gamm ./ a.cfd[idx].gamm[1]
    end
end

function drho(a::uwreal, mcid::Int64)
    idx = find_mcid(a, mcid)
    if (idx == 0)
        return [0.0]
    else
        return a.cfd[idx].drho
    end
end

function read_bdio(fb, ws::wspace)

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
        add_DB(dfl, convert(Int64, ids[i]), convert(Vector{Int64}, ivrep[is:ie]), ws)
        
        is = ie + 1
    end

    return uwreal(dfoo[1], p, d)
end

function write_bdio(p::uwreal, fb, iu::Int, ws::wspace)

    BDIO.BDIO_start_record!(fb, BDIO.BDIO_BIN_GENERIC, convert(Int32, iu), true)
    nid = convert(Int32, unique_ids!(p, ws))

    BDIO.BDIO_write!(fb, [p.mean], true)
    BDIO.BDIO_write!(fb, [nid], true)
    for j in 1:nid
        nd = convert(Int32, ws.fluc[ws.map_ids[p.ids[j]]].nd)
        BDIO.BDIO_write!(fb, [nd], true)
    end

    for j in 1:nid
        nr = convert(Int32, length(ws.fluc[ws.map_ids[p.ids[j]]].ivrep))
        BDIO.BDIO_write!(fb, [nr], true)
    end

    for j in 1:nid
        BDIO.BDIO_write!(fb, convert(Vector{Int32}, ws.fluc[ws.map_ids[p.ids[j]]].ivrep), true)
    end

    BDIO.BDIO_write!(fb, convert(Vector{Int32}, p.ids), true)
    
    for j in 1:nid
        for i in 1:length(p.prop)
            if (p.prop[i] && (ws.map_nob[i] == p.ids[j]))
                (nt,ip) = findmax(ws.fluc[i].ivrep)
                nt = convert(Int32, div(nt, 2*ws.fluc[i].ibn))
                BDIO.BDIO_write!(fb, [nt], true)
                continue
            end
        end
    end

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

    BDIO.BDIO_write_hash!(fb)

    return true
end

"""
    details(a::uwreal; io::IO=stdout, names::Dict{Int64, String} = Dict{Int64, String}())

Write out a detailed information on the error of `a`.

## Arguments

Optionally one can pass as a keyword argument (`io`) the `IO` stream to write to and a dictionary (`names`) that translates ensemble `ID`s into more human friendly `Strings`

## Example
```@example
using ADerrors # hide
a = uwreal(rand(2000),   120)
b = uwreal([1.2, 0.023], 121)
c = uwreal([5.2, 0.03],  122)
d = a + b - c
uwerr(d)

bnm = Dict{Int64, String}()
bnm[120] = "Very important ensemble"
bnm[122] = "H12B K87"
details(d, bnm)
```
"""
function details(a::uwreal, ws::wspace, io::IO=stdout, names::Dict{Int64, String} = Dict{Int64, String}())
    
    if (length(a.prop) == 0)
        print(a.mean)
        return
    end

    if (length(a.cfd) > 0) 
        println(a.mean, " +/- ", a.err)
        println(" ## Number of error sources: ", length(a.ids))

        n = 0
        for i in 1:length(a.cfd)
            if (length(a.cfd[i].gamm) > 0)
                n = n + 1
            end
        end
        println(" ## Number of MC ids       : ", n)
        println(" ## Contribution to error  :               Ensemble  [%]     [MC length]")

        
        truncate_ascii(s::String,n::Int) = s[1:min(sizeof(s),n)]
        ntrunc = 45
        v = zeros(length(a.cfd))
        for i in eachindex(v)
            v[i] = a.cfd[i].var
        end
        ip = sortperm(v, rev=true)
        for i in 1:length(a.cfd)
            sid = truncate_ascii(get(names, a.ids[ip[i]], string(a.ids[ip[i]])), ntrunc)
            if (length(a.cfd[ip[i]].gamm) > 0)
                idx  = ws.map_ids[a.ids[i]]
                nd   = ws.fluc[idx].nd
                Printf.@printf("  #  %45s %6.2f   %10d\n",
                        sid, 100.0 .* a.cfd[ip[i]].var ./ a.err^2, nd)
            else
                Printf.@printf("  #  %45s %6.2f            -\n",
                        sid, 100.0 .* a.cfd[ip[i]].var ./ a.err^2)
            end
        end
    else
        print(a.mean, " (Error not available... maybe run uwerr)")
    end
end

details(a::uwreal; io::IO=stdout, names::Dict{Int64, String} = Dict{Int64, String}()) = details(a, wsg, io, names)

read_uwreal(fb)  = read_bdio(fb, ADerrors.wsg)
write_uwreal(p::uwreal, fb, iu::Int) = write_bdio(p::uwreal, fb, iu::Int, ADerrors.wsg)

