###
### "THE BEER-WARE LICENSE":
### Alberto Ramos wrote this file. As long as you retain this 
### notice you can do whatever you want with this stuff. If we meet some 
### day, and you think this stuff is worth it, you can buy me a beer in 
### return. <alberto.ramos@cern.ch>
###
### file:    ADerrorsCF.jl
### created: Wed Jun 17 13:19:13 2020
###                               

const MIN_LENGTH = 500
const MIN_REP_LENGTH = 4
const nprm=50
const iprm = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
              43, 47, 53, 59, 61, 67, 71, 73,79,83,89,97,101,
              103,107,109,113,127,131,137,139,149,151,157,163,
              167,173,179,181,191,193,197,199,211,223,227,229]
const DEFAULT_STAU = 4.0
const DO_BIN       = false
const MINW         = 6

is_int32(str::String) =  Base.tryparse(Int32, str) !== nothing

function get_nbin_vec(nd::Array{Int64, 1}, dobin=DO_BIN)::Int64

    nbin::Int64 = 1
    if (!dobin)
        return nbin
    end
    if (any(nd .== 1))
        return nbin
    end
    
    nl = deepcopy(nd)
    for i = 1:nprm
        while (all(nl .% iprm[i] .== 0))
            nl = map(x->div(x,iprm[i]), nl)
            
            if (all(nl .< MIN_LENGTH))
                return nbin
            end
            nbin = nbin * iprm[i]
        end
    end
    return nbin
end

function bin_data(vin::Array{Float64, 1}, nbin::Int64)

    vout = similar(Array{Float64, 1}, div(length(vin), nbin))
    is::Int64 = 1
    for i = 1:length(vout)
        ie = is + nbin - 1
        vout[i] = Statistics.mean(vin[is:ie])
        is = ie+1
    end
    
    return vout
end

get_new_id(ws::wspace) = ws.newid = ws.newid - 1

function get_id_from_name(str::String, ws::wspace)

    if (haskey(ws.str2id, str))
        id = ws.str2id[str]
    else
        if is_int32(str)
            id = convert(Int64, Base.tryparse(Int32, str))
        else
            id = get_new_id(ws)
        end
        ws.id2str[id]  = str
        ws.str2id[str] = id
    end

    return id
end

function get_name_from_id(id::Int64, ws::wspace)

    if haskey(ws.id2str, id)
        str = ws.id2str[id]
    else
        str = string(id)
        ws.id2str[id]  = str
        ws.str2id[str] = id
    end

    return str
end

function add_repnames(id::Int64, ws::wspace, rname::Vector{String},
                      ridc::Vector{Int64})

    if haskey(ws.repnam, id)
        for i in 1:length(rname)
            if rname[i] != ws.repnam[id][i]
                error("Mistmatch in replica names for ensemble: "*get_name_from_id(id,ws))
            end
            if ridc[i] != ws.repidc[id][i]
                error("Mistmatch in replica configuration index for ensemble: "
                      *get_name_from_id(id,ws)*" replicum: "*rname[i])
            end
        end
    else
        ws.repnam[id] = rname
        ws.repidc[id] = ridc
    end
        
end
    
get_repnames_from_id(id::Int64, ws::wspace) = ws.repnam[id]
get_repidc_from_id(id::Int64, ws::wspace)   = ws.repidc[id]

function add_maps(id::Int64, ws::wspace, iv::Vector{Int64})

    ws.nob += 1
    push!(ws.map_nob, id)
    if (!haskey(ws.map_ids, id))
        ws.map_ids[id] = ws.nob
    else
        if (sum(iv) != ws.fluc[ws.map_ids[id]].nd)
            println(stderr, "ID:         ", ws.id2str[id])
            println(stderr, "DB  length: ", ws.fluc[ws.map_ids[id]].nd)
            println(stderr, "obs length: ", sum(iv))
            error("Mistmatch in data length for the same ensemble ID")
        end
        if (iv != ws.fluc[ws.map_ids[id]].ivrep)
            println(stderr, "ID:          ", ws.id2str[id])
            println(stderr, "DB  replica: ", ws.fluc[ws.map_ids[id]].ivrep)
            println(stderr, "obs replica: ", iv)
            error("Mistmatch in replica vector for the same ensemble ID")
        end
    end

    return nothing
end

function add_DB(delta::Vector{Float64}, id::Int64, iv::Vector{Int64}, ws::wspace, do_maps::Bool = true)


    nd = length(delta)

    if (nd == 1)
        new  = fbd(1, 1, delta, [1], Dict{Int64,Vector{Complex{Float64}}}())
        push!(ws.fluc, new)
    else
        if (sum(iv) != length(delta))
            println(stderr, "ID:          ", id)
            ArgumentError("Sum of replica length does not match number of measurements")
        end

        fseries = Dict{Int64,Vector{Complex{Float64}}}()
        nbin = get_nbin_vec(iv)
        is = 1
        for i in 1:length(iv)
            ie = is + iv[i] - 1
            nbdt = div(iv[i], nbin)
            datapad = [bin_data(delta[is:ie], nbin);
                       zeros(Float64, 600*(div(2*nbdt+1, 600)+1)-nbdt) ]
            fseries[i] = FFTW.fft(datapad)
            is = ie + 1
        end

        new = fbd(length(delta), nbin,
                  delta,
                  iv,
                  fseries)
        push!(ws.fluc, new)
    end

    if do_maps
        if (nd == 1)
            add_maps(id, ws, [1])
        else
            add_maps(id, ws, iv)
        end
    end

    return nothing
end

function uwcls(data::Vector{Float64}, id::Int64, ws::wspace, iv::Vector{Int64})
    if (length(data) == 2)
        add_DB([data[2]], id, Vector{Int64}(), ws)
        
        p = [false for n in 1:ws.nob]
        p[end] = true
        d = [0.0 for n in 1:ws.nob]
        d[end] = 1.0
        return uwreal(data[1], 0.0, 0.0,
                      p, d, Vector{Int64}(), Vector{cfdata}())
    else
        avg = Statistics.mean(data)
        add_DB(data .- avg, id, iv, ws)

        p = [false for n in 1:ws.nob]
        p[end] = true
        d = [0.0 for n in 1:ws.nob]
        d[end] = 1.0
        return uwreal(avg, 0.0, 0.0,
                      p, d, Vector{Int64}(), Vector{cfdata}())
        
    end
end

function uwcls_gaps(data::Vector{Float64},
                    id::Int64, ws::wspace,
                    iv::Vector{Int64},
                    idm::Vector{Int64},
                    nms::Int64)
    if (nms < 4)
        ArgumentError("MC data length has to be larger than 4")
    end
    if (sum(iv) != nms)
        ArgumentError("Sum of replica length does not match number of measurements")
    end

    avg = Statistics.mean(data)
    dt  = zeros(Float64, nms)
    for n in 1:length(idm)
        dt[idm[n]] = data[n] - avg
    end
    dt .= (nms/length(idm)) .* dt
    
    add_DB(dt, id, iv, ws)

    p = [false for n in 1:ws.nob]
    p[end] = true
    d = [0.0 for n in 1:ws.nob]
    d[end] = 1.0
    return uwreal(avg, 0.0, 0.0,
                  p, d, Vector{Int64}(), Vector{cfdata}())
end

function unique_ids!(a::uwreal, ws::wspace)

    if (length(a.ids) == 0)
        for i in 1:length(a.prop)
            if (a.prop[i])
                if (any(a.ids[1:end] .== ws.map_nob[i]))
                    continue
                end
                push!(a.ids, ws.map_nob[i])
            end
        end
    end
    return length(a.ids)
    
end

function wopt_ulli(nd::Int64, stau::Float64, gmm::Vector{Float64})

    tiw = 0.5
    if (gmm[1] == 0.0)
        return 1
    else
        for i in 2:length(gmm)
            tiw = tiw + gmm[i]/gmm[1]
            if (tiw <= 0.5)
                return max(MINW,i)
            else
                tau = stau/log((2.0*tiw+1.0)/(2.0*tiw-1.0))
                gw = exp(-(i-1.0)/tau) - tau/sqrt((i-1.0)*nd)
                if (gw < 0.0)
                    return max(MINW,i)
                end
            end
        end
    end
    
end

function uwerror(a::uwreal, ws::wspace, wpm::Dict{Int64,Vector{Float64}})

    nid = unique_ids!(a, ws)

    if (length(a.cfd) == 0)
        a.cfd = Vector{cfdata}(undef, nid)
        for j in 1:nid
            a.cfd[j] = cfdata()
        end
    end
    
    for j in 1:nid
        wp = get(wpm, a.ids[j], [-1.0,-1.0,-1.0,-1.0,1.0])
        if length(wp) == 4
            nbin = 1
        else
            nbin = round(Int64, wp[5])
            if nbin < 1
                nbin = 1
            end
        end
        if (a.cfd[j].iw != 0) && (nbin == a.cfd[j].nbin)
            continue
        end

        a.cfd[j].nbin = nbin
        idx = ws.map_ids[a.ids[j]]
        nd  = ws.fluc[idx].nd
        if (nd != 1)
            nd_eff = div(nd, nbin)
            (nt,ip) = findmax(ws.fluc[idx].ivrep)
            nt = div(nt, 2*a.cfd[j].nbin)
            nrep   = length(ws.fluc[idx].ivrep)
            ftemp  = Dict{Int64,Vector{Complex{Float64}}}()
            @inbounds for i in 1:nrep
                ns = length(ws.fluc[idx].fourier[i])
                ftemp[i] = zeros(Complex{Float64}, ns)
            end
        end
        
        @inbounds for i in 1:length(a.prop)
            if (a.prop[i] && (ws.map_nob[i] == a.ids[j]))
                if (nd == 1)
                    a.cfd[j].var = a.cfd[j].var + a.der[i]*ws.fluc[i].delta[1]
                else
                    @inbounds for k in 1:nrep # ensemble growth here
                        ftemp[k] = ftemp[k] + a.der[i]*ws.fluc[i].fourier[k]
                    end
                end
            end
        end
        
        if (nd == 1)
            a.cfd[j].var = a.cfd[j].var^2
            a.cfd[j].iw  = 1
        else
            a.cfd[j].gamm = zeros(nt)
            @inbounds for k in 1:nrep # binning here
                if nbin == 1
                    ftemp[k] .= ftemp[k].*conj(ftemp[k])
                    FFTW.ifft!(ftemp[k])
                    
                    @inbounds for ig in 1:min(nt,length(ftemp[k]))
                        a.cfd[j].gamm[ig] = a.cfd[j].gamm[ig] + real(ftemp[k][ig])
                    end
                else
                    println(ws.fluc[idx].ivrep[k])
                    ns = length(ftemp[k])
                    if ns % nbin != 0
                        if nd % nbin != 0
                            error("Bin size ($nbin) not allowed for replica length (",
                                  ws.fluc[idx].ivrep[k], ")")
                        end
                        ftt = FFTW.ifft(ftemp[k])
                        ns = nbin*(div(2*nd+1, nbin)+1)
                        resize!(ftt, ns)
                        @inbounds for kk in length(ftemp[k])+1:ns
                            ftt[kk] = zero(eltype(ftt))
                        end
                        FFTW.fft!(ftt)
                        nse = div(ns,nbin)
                        
                        ft2 = zeros(eltype(ftt), nse)
                        @inbounds for kk in 2:nse
	                    @inbounds for r in 1:nbin
	                        zr = exp(im  * 2*pi*(kk + (r-1)*nse - 1)/ns)
	                        ft2[kk] += (ftt[kk+(r-1)*nse]/nbin^2) * (1-zr^nbin)/(1-zr)
	                    end
                        end
                    else
                        nse = div(ns,nbin)
                        
                        ft2 = zeros(eltype(ftemp[k]), nse)
                        @inbounds for kk in 2:nse
	                    @inbounds for r in 1:nbin
	                        zr = exp(im  * 2*pi*(kk + (r-1)*nse - 1)/ns)
	                        ft2[kk] += (ftemp[k][kk+(r-1)*nse]/nbin^2) * (1-zr^nbin)/(1-zr)
	                    end
                        end
                    end

                    ft2 .= ft2 .* conj.(ft2)
                    FFTW.ifft!(ft2)
                    @inbounds for ig in 1:min(nt,length(ft2))
                        a.cfd[j].gamm[ig] = a.cfd[j].gamm[ig] + real(ft2[ig])
                    end
                end
            end
            @inbounds for ig in 1:nt
                nrcnt = 0
                for k in 1:nrep
                    if (div(ws.fluc[idx].ivrep[k], 2*ws.fluc[idx].ibn) > ig-1 )
                        nrcnt = nrcnt + 1
                    end
                end
                a.cfd[j].gamm[ig] = a.cfd[j].gamm[ig] / (nd_eff - nrcnt*(ig-1))
            end
            
            iw = wopt_ulli(nd_eff, DEFAULT_STAU, a.cfd[j].gamm)
            a.cfd[j].iw = iw
            
            gwin  = view(a.cfd[j].gamm, 2:iw)
#            dbias = a.cfd[j].gamm[1] + 2.0*sum(gwin)
#            if (dbias > 0.0)
#                a.cfd[j].gamm .= a.cfd[j].gamm .+ dbias/nd_eff
#            end
            
            a.cfd[j].drho = zeros(nt)
            if (a.cfd[j].gamm[1] != 0.0)
                @inbounds for i in 1:nt
                    is = max(1, i-iw-2) + 1
                    ie = is + iw-1
                    @inbounds for k in is:ie
                        if (k < nt+1)
                            cont = -2.0*a.cfd[j].gamm[k]*a.cfd[j].gamm[i]/a.cfd[j].gamm[1]^2
                        else
                            cont = 0.0
                        end
                        
                        if ((i+k-2) < nt)
                            cont = cont + a.cfd[j].gamm[i+k-1]/a.cfd[j].gamm[1]
                        end
                        if (abs(i-k) < nt)
                            cont = cont + a.cfd[j].gamm[abs(i-k)+1]/a.cfd[j].gamm[1]
                        end
                        a.cfd[j].drho[i] = a.cfd[j].drho[i] + cont^2
                    end
                    a.cfd[j].drho[i] = sqrt(a.cfd[j].drho[i]/nd_eff)
                end
            else
                a.cfd[j].var = a.cfd[j].var^2
            end
        end
    end

    id_neg_taui = Vector{Int64}()
    a.err  = 0.0
    a.derr = 0.0
    for j in 1:nid
        idx = ws.map_ids[a.ids[j]]
        if (ws.fluc[idx].nd == 1)
            a.cfd[j].taui  = 0.5
            a.cfd[j].dtaui = 0.0
            vti = 0.0
            iw  = 0
        else
            nd_eff = div(ws.fluc[idx].nd, a.cfd[j].nbin)
            (nt,ip) = findmax(ws.fluc[idx].ivrep)
            nt = div(nt, 2*a.cfd[j].nbin)
            if (a.cfd[j].gamm[1] == 0.0)
                a.cfd[j].taui  = 0.0
                a.cfd[j].dtaui = 0.0
                vti = 0.0
                iw  = 0
            else
                wp = zeros(4)
                if haskey(wpm, a.ids[j]) 
                    wp = get(wpm, a.ids[j], [-1.0,-1.0,-1.0,-1.0,1.0])
                    if (wp[1] > 0.0)
                        a.cfd[j].iw = round(Int64, wp[1])
                    elseif (wp[2] > 0.0)
                        a.cfd[j].iw = wopt_ulli(nd_eff, wp[2], a.cfd[j].gamm)
                    end
                    
                    if (wp[3] > 0.0)
                        iw = 1
                        for k in 2:nt
                            if (a.cfd[j].drho[k]*wp[3] > a.cfd[j].gamm[k]/a.cfd[j].gamm[1])
                                iw = k-1
                                break
                            end
                        end
                        a.cfd[j].iw = iw
                    end
                    
                    if (wp[4] > 0.0)
                        texp = wp[4]/a.cfd[j].nbin
                    else
                        texp = 0.0
                    end
                else
                    texp = 0.0
                    a.cfd[j].iw = wopt_ulli(nd_eff, DEFAULT_STAU, a.cfd[j].gamm)
                end
                iw = a.cfd[j].iw
                    
                gwin  = view(a.cfd[j].gamm, 2:iw)
                vti = 0.5 + sum(gwin)/a.cfd[j].gamm[1]
                a.cfd[j].dtaui = sqrt(vti^2 * (4.0*iw-2.0*vti+2.0)/nd_eff)
                a.cfd[j].taui  = vti + texp*a.cfd[j].gamm[iw+1]/a.cfd[j].gamm[1]
            end
            
            if (a.cfd[j].taui >= 0.0)
                a.cfd[j].var = a.cfd[j].gamm[1] * 2.0*a.cfd[j].taui/nd_eff
            else
                push!(id_neg_taui, a.ids[j])
            end
        end

        a.err = a.err + a.cfd[j].var
        if (iw > 1)
            if vti*(a.cfd[j].iw-0.5)/nd_eff > 0.0
                a.derr = a.derr + vti*(a.cfd[j].iw-0.5)/nd_eff
            end
        end
    end

    if (length(id_neg_taui) == 0)
        a.err  = sqrt(a.err)
        a.derr = sqrt(a.derr)
    else
        println(a.cfd[1].gamm[1:a.cfd[1].iw+1])
        println(id_neg_taui)
        println(stderr, "ID's with negative tau_int: ")
        for i in 1:length(id_neg_taui)
            println("     ",  get_name_from_id(id_neg_taui[i], ws))
        end
        error("Error analysis failed for some ID's. Choose your window more carefully")
    end

    return nothing
end

function unique_ids_multi(a::Vector{uwreal}, ws::wspace)

    is = 0
    for i in 1:length(a)
        is = is + unique_ids!(a[i], ws)
    end

    n = 0
    ids = Vector{Int64}()
    for k in 1:length(a)
        for i in 1:length(a[k].prop)
            if (a[k].prop[i])
                if (any(ids[1:end] .== ws.map_nob[i]))
                    continue
                end
                n = n + 1
                push!(ids, ws.map_nob[i])
            end
        end
    end

    return ids
end


function cov(a::Vector{uwreal}, ws::wspace, wpm::Dict{Int64,Vector{Float64}})

    for i in 1:length(a)
        ADerrors.uwerror(a[i], ws, wpm)
    end
    ids = unique_ids_multi(a, ws)
    nid = length(ids)
    iw  = zeros(Int64, nid)

    for k in 1:length(a)
        for j in 1:length(a[k].ids)
            for i in 1:nid
                if (a[k].ids[j] == ids[i])
                    iw[i] = max(iw[i], a[k].cfd[j].iw)
                end
            end
        end
    end
    
    wopt = Dict{Int64,Vector{Float64}}()
    for i in 1:nid
        wopt[ids[i]] = [Base.convert(Float64, iw[i]), -1.0, -1.0, -1.0]
    end
    
    cov = zeros(Float64, length(a), length(a))
    for k in 1:length(a)
        ADerrors.uwerror(a[k], ws, wopt)
        cov[k,k] = a[k].err^2
    end
    
    for j in 1:nid
        idx  = ws.map_ids[ids[j]]
        nd   = ws.fluc[idx].nd
        if (nd != 1)
            nd_eff = div(nd, ws.fluc[idx].ibn)
            (nt,ip) = findmax(ws.fluc[idx].ivrep)
            nt = div(nt, 2*ws.fluc[idx].ibn)
            nrep = length(ws.fluc[idx].ivrep)
            
            ftemp1 = Dict{Int64,Vector{Complex{Float64}}}()
            ftemp2 = Dict{Int64,Vector{Complex{Float64}}}()
            for i in 1:nrep
                ns = length(ws.fluc[idx].fourier[i])
                ftemp1[i] = zeros(Complex{Float64}, ns)
                ftemp2[i] = zeros(Complex{Float64}, ns)
            end

            gamm = zeros(Float64, nt)
        end
        for k1 in 1:length(a)-1
            v1 = 0.0
            for i in 1:length(a[k1].prop)
                if (a[k1].prop[i] && (ws.map_nob[i] == ids[j]))
                    if (nd == 1)
                        v1 = v1 + a[k1].der[i]*ws.fluc[i].delta[1]
                    else
                        for k in 1:nrep
                            ftemp1[k] = ftemp1[k] + a[k1].der[i]*ws.fluc[i].fourier[k]
                        end
                    end
                end
            end
            for k2 in k1+1:length(a)
                v2 = 0.0
                for i in 1:length(a[k2].prop)
                    if (a[k2].prop[i] && (ws.map_nob[i] == ids[j]))
                        if (nd == 1)
                            v2 = v2 + a[k2].der[i]*ws.fluc[i].delta[1]
                        else
                            for k in 1:nrep
                                ftemp2[k] = ftemp2[k] + a[k2].der[i]*ws.fluc[i].fourier[k]
                            end
                        end
                    end
                end

                if (nd == 1)
                    cov[k1,k2] = cov[k1,k2] + v1*v2
                else
                    for k in 1:nrep
                        ftemp2[k] .= ftemp1[k].*conj(ftemp2[k])
                        FFTW.ifft!(ftemp2[k])
                        
                        for ig in 1:min(nt,length(ftemp2[k]))
                            gamm[ig] = gamm[ig] + real(ftemp2[k][ig])
                        end
                    end

                    for ig in 1:nt
                        nrcnt = count(map(x -> div(x, 2*ws.fluc[idx].ibn), ws.fluc[idx].ivrep) .> ig-1)
                        gamm[ig] = gamm[ig] / (nd_eff - nrcnt*(ig-1))
                    end
                    
                    dbias = gamm[1] + 2.0*sum(gamm[2:iw[j]])
                    gamm .= gamm .+ dbias/nd_eff

                    cov[k1,k2] = cov[k1, k2] + ( gamm[1] + 2.0*sum(gamm[2:iw[j]]) )/nd_eff
                
                    for i in 1:nrep
                        for l in 1:length(ws.fluc[idx].fourier[i])
                            ftemp2[i][l] = Base.convert(Complex{Float64}, 0.0)
                        end
                    end
                    for i in 1:nt
                        gamm[i] = 0.0
                    end
                end
            end
            
            if (nd > 1)
                for i in 1:nrep
                    for l in 1:length(ws.fluc[idx].fourier[i])
                        ftemp1[i][l] = Base.convert(Complex{Float64}, 0.0)
                    end
                end
            end
        end
    end

    for k1 in 1:length(a)-1
        for k2 in k1+1:length(a)
            cov[k2,k1] = cov[k1,k2]
        end
    end
    
    return cov
end

function trcov(M, a::Vector{uwreal}, ws::wspace, wpm::Dict{Int64,Vector{Float64}})

    usvt = LinearAlgebra.svd(M)
    n    = length(a)

    p = similar(a)
    for i in 1:n
        p[i] = a[1]*usvt.U[1,i]
        for j in 2:n
            p[i] = p[i] + a[j]*usvt.U[j,i]
        end
        ADerrors.uwerror(p[i], ws, wpm)
    end

    ids = unique_ids_multi(a, ws)
    iw = zeros(Int64, length(ids))
    for k in 1:n
        for j in 1:length(p[k].ids)
            for i in 1:length(ids)
                if (p[k].ids[j] == ids[i])
                    iw[i] = max(iw[i], p[k].cfd[j].iw)
                end
            end
        end
    end
    
    wopt = Dict{Int64,Vector{Float64}}()
    for i in 1:length(ids)
        wopt[ids[i]] = [Base.convert(Float64, iw[i]), -1.0, -1.0, -1.0]
    end

    tr = 0.0
    for k in 1:n
        ADerrors.uwerror(p[k], ws, wopt)
        tr = tr + usvt.S[k]*p[k].err^2
    end

    return tr
end

function trcorr(M, a::Vector{uwreal}, ws::wspace, wpm::Dict{Int64,Vector{Float64}}, W::Vector{Float64})

    usvt = LinearAlgebra.svd(M)
    n    = length(a)

    Ww = zeros(Float64, length(a))
    if (length(W) == 0)
        for i in 1:length(a)
            ADerrors.uwerror(a[i], ws, wpm)
            Ww[i] = a[i].err
        end
    else
        Ww .= W
    end
    
    p = similar(a)
    for i in 1:n
        p[i] = a[1]*usvt.U[1,i]/Ww[1]
        for j in 2:n
            p[i] = p[i] + a[j]*usvt.U[j,i]/Ww[j]
        end
        ADerrors.uwerror(p[i], ws, wpm)
    end

    ids = unique_ids_multi(a, ws)
    iw = zeros(Int64, length(ids))
    for k in 1:n
        for j in 1:length(p[k].ids)
            for i in 1:length(ids)
                if (p[k].ids[j] == ids[i])
                    iw[i] = max(iw[i], p[k].cfd[j].iw)
                end
            end
        end
    end
    
    wopt = Dict{Int64,Vector{Float64}}()
    for i in 1:length(ids)
        wopt[ids[i]] = [Base.convert(Float64, iw[i]), -1.0, -1.0, -1.0]
    end

    tr = 0.0
    for k in 1:n
        ADerrors.uwerror(p[k], ws, wopt)
        tr = tr + usvt.S[k]*p[k].err^2
    end

    return tr
end

##                                         ##
# Module workspace and function definitions #
##                                         ##


wsg = ADerrors.wspace(similar(Vector{ADerrors.fbd}, 0),
                      0,
                      similar(Vector{Int64}, 0),
                      Dict{Int64, Int64}(),
                      Dict{Int64, String}(), Dict{String, Int64}(),
                      Dict{Int64, Vector{String}}(),
                      Dict{Int64, Vector{Int64}}(),
                      -12345)

get_id_from_name(str::String) = get_id_from_name(str, wsg)
taui(a::uwreal, str::String) = taui(a, get_id_from_name(str))
dtaui(a::uwreal, str::String) = dtaui(a, get_id_from_name(str))
window(a::uwreal, str::String) = window(a, get_id_from_name(str))
rho(a::uwreal, str::String) = rho(a, get_id_from_name(str))
drho(a::uwreal, str::String) = drho(a, get_id_from_name(str))
mchist(a::uwreal, str::String) = mchist(a, get_id_from_name(str))
err(a::uwreal, str::String) = err(a, get_id_from_name(str))

empt = Dict{Int64,Vector{Float64}}()

function dict_name_to_id(wpm::Dict{String, Vector{Float64}})
    wp = Dict{Int64, Vector{Float64}}()
    for i in keys(wpm)
        wp[get_id_from_name(i)] = wpm[i]
    end
    
    return wp
end

"""
    change_id(; from::String, to::String)

Changes the ensemble `id` from the string/integer `from` to the string `to`. Note that a call to this routine will affect the output of all currently defined `uwreal` types that depend on ensemble `ID` `from`. On the other hand, it will not affect data that is read after the call to `change_id`. This means that:
- All input should be performed **before** any calls to `change_id`.
- All output (i.e. calls to `write_uwreal` or `details`) is affected by previous calls to `change_id`

```@example
using ADerrors # hide
a = uwreal([1.2, 0.2], "Simple var with error")   # a = 1.2 +/- 0.2
b = uwreal([1.2, 0.2], 1233)  # c = 1.2 +/- 0.2

d = a-b
uwerr(d)
details(d)

change_id(from=1233, to="Error source from experiment")
change_id(from="Simple var with error", to="Error source from simulations")
details(d)
```
"""
function change_id(ws::wspace; from::Union{String,Int64}, to::String)

    if (isa(from, Int64))
        id = from
    else
        id = get(ws.str2id, from, nothing)
    end
    
    if id == nothing
        error("ID does not exists in database")
    else
        if haskey(ws.str2id, to)
            error("Destination ID already exists in the database")
        else
            ws.id2str[id] = to
            ws.str2id[to] = id
            
            if (isa(from, Int64))
                delete!(ws.str2id, string(from))
            else
                delete!(ws.str2id, from)
            end
        end
    end

    return nothing
end

change_id(; from::Union{String,Int64}, to::String) = change_id(wsg, from=from, to=to)


"""
    uwreal(x::Float64)
    uwreal([value::Float64, error::Float64], mcid)
    uwreal(data::Vector{Float64}, mcid[, replica::Vector{Int64}])
    uwreal(data::Vector{Float64}, mcid[, replica::Vector{Int64}], idm::Vector{Int64}, nms::Int64)

Returns an `uwreal` data type. Depending on the first argument, the `uwreal` stores the following information:
- Input is a single `Float64`. In this case the variable acts in exactly the same way as a real number. This is understood as a quantity with zero error.
- Input is a 2 element Vector of `Float64` `[value, error]`. In this case the data is understood as `value +/- error`.
- Input is a Vector of `Float64` of length larger than 4. In this case the data is understood as consecutive measurements of an observable in a Monte Carlo (MC) simulation. 

In the last two cases, an ensemble `ID` is required as input. Data with the same `ID` are considered as correlated (i.e. fully correlated for the case of a `value +/- error` observables measured on the same sample for the case of data from a MC simulation). The preferred way to input the ensemble tag is via a `String` that uniquely identifies the ensemble, but an integer is also supported for legacy reasons. For example:

```@example
using ADerrors # hide
a = uwreal([1.2, 0.2], "Simple var with error")   # a = 1.2 +/- 0.2
b = uwreal([1.2, 0.2], "Simple var with error")   # b = 1.2 +/- 0.2
c = uwreal([1.2, 0.2], "Another var with error")  # c = 1.2 +/- 0.2

d = a-b
uwerr(d)
println("d has zero error because a and b are correlated", d)

e = a-c
uwerr(e)
println("e has non zero error because a and c are independent", e)
```

### Replica

`data` can contain measurements in several replica (i.e. independent simulations with the same physical and algorithmic parameters). How many measurements correspond to each replica are specified by an optional replica vector argument `replica`. Note that `length(replica)` is the number of replica (i.e. independent simulations), and that `sum(replica)` must match `length(data)` 
```@example
using ADerrors # hide
# 1000 measurements in three replica of lengths
# 500, 100 and 400
a = uwreal(rand(1000), "Ensemble with three replica", [500, 100, 400]) 
```
### Gaps in the measurements

In some situations an observable is not measured in every configuration. In this case two additional arguments are needed to define the observable
- `idm`. Type `Vector{Int64}`. `idm[n]` labels the configuration where `data[n]`is measured.
- `nms`. Type `Int64`. The total number of measurements in the ensemble
```@example
using ADerrors # hide
# Observable measured on the odd configurations 
# 1, 3, 5, ..., 999 on an emsemble of length 1000
a = uwreal(rand(500), "Observable with gaps", collect(1:2:999), 1000)

# Observable measured on the first 900 configurations 
# on the same emsemble
b = uwreal(rand(900), "Observable with gaps", collect(1:900), 1000)
```
Note that in this case, if the ensemble has different replica, `sum(replica)` must match `nms`
```@example
using ADerrors # hide
# Observable measured on the even configurations 
# 2, 4, 6, ..., 200 on an emsemble of length 200
# with two replica of lengths 75, 125
a = uwreal(data_a[1:500], "Observable with gaps in an ensemble with replica", [75, 125], collect(2:2:200), 200)
```
"""
uwreal(x::Float64) = ADerrors.uwreal(x, 0.0, 0.0, 
                                     Vector{Bool}(), Vector{Float64}(), 
                                     Vector{Int64}(), Vector{cfdata}())

uwreal(data::Vector{Float64},
       id::Int64) = uwreal(data,
                           string(id))
uwreal(data::Vector{Float64},
       id::Int64,
       iv::Vector{Int64}) = uwreal(data,
                                   string(id),
                                   iv)
uwreal(data::Vector{Float64},
       id::Int64,
       idm::Vector{Int64},
       nms::Int64) = uwreal(data,
                            string(id),
                            [nms],
                            idm,
                            nms)
uwreal(data::Vector{Float64},
       id::Int64,
       iv::Vector{Int64},
       idm::Vector{Int64},
       nms::Int64) = uwreal(data,
                            string(id),
                            iv,
                            idm,
                            nms)

function uwreal(data::Vector{Float64}, str::String)
    uw = ADerrors.uwcls(data,
                        get_id_from_name(str, wsg),
                        wsg,
                        [length(data)])
    v = [str*"_r0"]
    idc = collect(1:length(data))
    add_repnames(get_id_from_name(str, wsg), wsg, v, idc)

    return uw
end

function uwreal(data::Vector{Float64},
                str::String,
                iv::Vector{Int64})
    
    uw = ADerrors.uwcls(data, get_id_from_name(str, wsg), wsg, iv)
    v = Vector{String}(undef, length(iv))
    idc = Vector{Int64}(undef, length(data))
    iof = 0
    for i in 1:length(v)
        v[i] = str*"_r"*string(i-1)
        for k in 1:iv[i]
            idc[k+iof] = k
        end
        iof = iof + iv[i]
    end

    add_repnames(get_id_from_name(str, wsg), wsg, v, idc)
    
    return uw
end

function uwreal(data::Vector{Float64},
                str::String,
                idm::Vector{Int64},
                nms::Int64)

    uw = ADerrors.uwcls_gaps(data,
                    get_id_from_name(str, wsg), wsg,
                    [nms],
                    idm,
                    nms)
    v = [str*"_r0"]
    idc = collect(1:nms)
    add_repnames(get_id_from_name(str, wsg), wsg, v, idc)

    return uw
end

function uwreal(data::Vector{Float64},
                str::String,
                iv::Vector{Int64},
                idm::Vector{Int64},
                nms::Int64)
    uw = ADerrors.uwcls_gaps(data,
                             get_id_from_name(str, wsg), wsg, 
                             iv,
                             idm,
                             nms)
    v = Vector{String}(undef, length(iv))
    idc = Vector{Int64}(undef, nms)
    iof = 0
    for i in 1:length(v)
        v[i] = str*"_r"*string(i-1)
        for k in 1:iv[i]
            idc[k+iof] = k
        end
        iof = iof + iv[i]
    end
    add_repnames(get_id_from_name(str, wsg), wsg, v, idc)

    return uw
end
    

function uwreal(data::Vector{Float64}, str::String, rname::Vector{String})
    uw = ADerrors.uwcls(data,
                        get_id_from_name(str, wsg),
                        wsg,
                        [length(data)])
    idc = collect(1:length(data))
    add_repnames(get_id_from_name(str, wsg), wsg, rname, idc)

    return uw
end

function uwreal(data::Vector{Float64},
                str::String, rname::Vector{String},
                iv::Vector{Int64})
    
    uw = ADerrors.uwcls(data, get_id_from_name(str, wsg), wsg, iv)

    idc = Vector{Int64}(undef, length(data))
    iof = 0
    for i in 1:length(iv)
        for k in 1:iv[i]
            idc[k+iof] = k
        end
        iof = iof + iv[i]
    end
    add_repnames(get_id_from_name(str, wsg), wsg, rname, idc)
    
    return uw
end

function uwreal(data::Vector{Float64},
                str::String, rname::Vector{String},
                idm::Vector{Int64},
                nms::Int64)

    uw = ADerrors.uwcls_gaps(data,
                    get_id_from_name(str, wsg), wsg,
                    [nms],
                    idm,
                    nms)

    idc = collect(1:nms)
    add_repnames(get_id_from_name(str, wsg), wsg, rname, idc)

    return uw
end

function uwreal(data::Vector{Float64},
                str::String, rname::Vector{String},
                iv::Vector{Int64},
                idm::Vector{Int64},
                nms::Int64)
    uw = ADerrors.uwcls_gaps(data,
                             get_id_from_name(str, wsg), wsg, 
                             iv,
                             idm,
                             nms)

    idc = Vector{Int64}(undef, nms)
    iof = 0
    for i in 1:length(iv)
        for k in 1:iv[i]
            idc[k+iof] = k
        end
        iof = iof + iv[i]
    end
    add_repnames(get_id_from_name(str, wsg), wsg, rname, iv)

    return uw
end


@doc raw"""
     uwerr(a::uwreal[, wpm::Dict{Int64,Vector{Float64}}])
     uwerr(a::uwreal[, wpm::Dict{String,Vector{Float64}}])

Performs error analysis on the observable `a`.
```@example
using ADerrors # hide
a = uwreal([1.3, 0.01], "Var with error") # 1.3 +/- 0.01
b = sin(2.0*a)
uwerr(b)
println("Look how I propagate errors:              ", b)
c = 1.0 + b - 2.0*sin(a)*cos(a)
uwerr(c)
println("Look how good I am at this (zero error!): ", c)
```

### Optimal window

Error in data coming from a Monte Carlo ensemble is determined by summing the autocorrelation function ``\Gamma_{\rm ID}(t)`` of the data for each ensemble ID. In practice this sum is truncated up to a window ``W_{\rm ID}``.

By default, the summation window is determined as [proposed by U. Wolff](https://inspirehep.net/literature/621085) with a parameter ``S_{\tau} = 4``, but other methods are avilable via the optional argument `wpm`. 

For each ensemble, one a pass the following parameters:
- `window`: The autocorrelation function is summed up to `t = round(vp[1])`.
- `stau`: The sumation window is determined [using U. Wolff poposal](https://inspirehep.net/literature/621085) with ``S_{\tau} = {\rm vp[2]}``.
- `signal`: The autocorrelation function ``\Gamma(t)`` is summed up a point where its error ``\delta\Gamma(t)`` is a factor `vp[3]` times larger than the signal. 

These three parameters fix the summation window. The following options also affect the error analysis
- `texp`: Add a tail to the autocorrelation function with this value of ``\tau_{\rm exp}``. See [the reference](https://inspirehep.net/literature/871175) for an explanation of the procedure.
- `bin`: Bin the data before computing the autocorrelation function, usg this value as bin size.
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

# Load the data in a uwreal
a = uwreal(x.^2, "Random walk in [-1,1]")
wpm = Dict{String,Dict{String, Real}}()

# Use default analysis (stau = 4.0)
uwerr(a)
println("default:                   ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")
println("  - window: ", window(a, "Random walk in [-1,1]"))

# This will still do default analysis because 
# a does not depend on emsemble foo
wpm["Ensemble foo"] = Dict("window" => 34)
uwerr(a, wpm)
println("default:                   ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")
println("  - window: ", window(a, "Random walk in [-1,1]"))

# Use default analysis (stau = 4.0), bin data with bin size 5
wpm["Random walk in [-1,1]"] = Dict("bin" => 5)
uwerr(a, wpm)
println("default (bin size 5):     ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")
println("  - window: ", window(a, "Random walk in [-1,1]"))

# Fix the summation window to 1 (i.e. uncorrelated analysis)
wpm["Random walk in [-1,1]"] = Dict("window" => 1)
uwerr(a, wpm)
println("uncorrelated:              ",  a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")
println("  - window: ", window(a, "Random walk in [-1,1]"))

# Use stau = 1.5
wpm["Random walk in [-1,1]"] = Dict("stau" => 1.5)
uwerr(a, wpm)
println("stau = 1.5:                ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")
println("  - window: ", window(a, "Random walk in [-1,1]"))

# Use fixed window 2 and add tail with texp = 100.0
wpm["Random walk in [-1,1]"] = Dict("window" => 2, "texp" => 100.0)
uwerr(a, wpm)
println("Fixed window 2, texp=100: ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")
println("  - window: ", window(a, "Random walk in [-1,1]"))

# Sum up to the point that the signal in Gamma is 
# 1.5 times the error and add a tail with texp = 10.0
wpm["Random walk in [-1,1]"] = Dict("signal" => 1.5, "texp" => 10.0)
uwerr(a, wpm)
println("signal/noise=1.5, texp=10: ", a, " (tauint = ", taui(a, "Random walk in [-1,1]"), ")")
println("  - window: ", window(a, "Random walk in [-1,1]"))
```


#### Old (legacy) way to set the window parameters

For each ensemble `ID` one can pass a vector of `Float64` of length 4. The first three components of the vector specify the criteria to determine the summation window:
- `vp[1]`: The autocorrelation function is summed up to `t = round(vp[1])`.
- `vp[2]`: The sumation window is determined [using U. Wolff poposal](https://inspirehep.net/literature/621085) with ``S_{\tau} = {\rm vp[2]}``.
- `vp[3]`: The autocorrelation function ``\Gamma(t)`` is summed up a point where its error ``\delta\Gamma(t)`` is a factor `vp[3]` times larger than the signal. 

An additional fourth parameter `vp[4]`, tells `ADerrors` to add a tail to the error with ``\tau_{\rm exp} = {\rm vp[4]}``. See [the reference](https://inspirehep.net/literature/871175) for an explanation of the procedure.

Note that:
- Negative values of `vp[1:4]` are ignored. 
- One, and only one, of the components `vp[1:3]` has to be positive. This chooses your criteria to determine the summation window.
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
```
"""
uwerr(a::uwreal, wpm::Dict{Int64,Vector{Float64}}) = ADerrors.uwerror(a, wsg, wpm)
uwerr(a::uwreal) = ADerrors.uwerror(a::uwreal, wsg, empt)
uwerr(a::uwreal, wpm::Dict{String,Vector{Float64}}) = uwerr(a, dict_name_to_id(wpm))
function uwerr(a::uwreal, wpm::Dict{String,Dict{String, Real}})
    
    wpi = Dict{Int64, Vector{Float64}}()
    for (k, v) in wpm
        r1 = convert(Float64, get(v, "window", -1))
        r2 = convert(Float64, get(v, "stau", -1))
        r3 = convert(Float64, get(v, "signal", -1))
        r4 = convert(Float64, get(v, "texp", -1))
        r5 = convert(Float64, get(v, "bin", 1))
        wpi[get_id_from_name(k)] = [r1,r2,r3,r4,r5]
    end
    
    return uwerr(a, wpi)
end


"""
    cov(a::Vector{uwreal}[, wpm])

Determine the covariance matrix between the vector of observables `a[:]`. 
```@example
using ADerrors, LinearAlgebra # hide
a = uwreal([1.3, 0.01], "Var 1") # 1.3 +/- 0.01
b = uwreal([5.3, 0.23], "Var 2") # 5.3 +/- 0.23
uwerr(a)
uwerr(b)

x = [a+b, a-b]
mat = ADerrors.cov(x)
println("Covariance: ", mat[1,1], " ", mat[1,2])
println("            ", mat[2,1], " ", mat[2,2])
println("Check (should be zero): ",  mat[1,1] - mat[2,2])
println("Check (should be zero): ",  mat[1,2] - (err(a)^2-err(b)^2))
```

# Case of Monte Carlo data

An optional parameter `wpm` can be used to choose the summation window for the relevant autocorrelation functions. The situation is completely analogous to the case of error analysis of single variables.
"""
cov(a::Vector{uwreal}) = cov(a::Vector{uwreal}, wsg, empt)
cov(a::Vector{uwreal}, wpm::Dict{Int64,Vector{Float64}}) = cov(a, wsg, wpm)
cov(a::Vector{uwreal}, wpm::Dict{String,Vector{Float64}}) = cov(a, wsg, dict_name_to_id(wpm))


@doc raw"""
     trcov(M::Array{Float64, 2}, a::Vector{uwreal}[, wmp])

Given a vector of `uwreal`, `a[:]` and a two dimensional symmtric positive definite array `M`, this routine computes  ``{\rm tr}(MC)``, where ``C_{ij} = {\rm cov}(a[i], a[j])``. 
```@example
using ADerrors, LinearAlgebra # hide
a = uwreal([1.3, 0.01], "Var with error 1") # 1.3 +/- 0.01
b = uwreal([5.3, 0.23], "Var with error 2") # 5.3 +/- 0.23
c = uwreal(rand(2000), "White noise ensemble")

x = [a+b+sin(c), a-b+cos(c), c-b/a]
M = [1.0 0.2 0.1
     0.2 2.0 0.3
     0.1 0.3 1.0]

mcov = cov(x)
d = tr(mcov * M)
println("Better be zero: ", d -trcov(M, x))

# Case of Monte Carlo data

An optional parameter `wpm` can be used to choose the summation window for the relevant autocorrelation functions. The situation is completely analogous to the case of error analysis of single variables.
```
"""
trcov(M, a::Vector{uwreal}) = trcov(M, a, wsg, empt)
trcov(M, a::Vector{uwreal}, wpm::Dict{Int64,Vector{Float64}}) = trcov(M, a, wsg, wpm)
trcov(M, a::Vector{uwreal}, wpm::Dict{String,Vector{Float64}}) =
    trcov(M, a, wsg, dict_name_to_id(wpm))

trcorr(M, a::Vector{uwreal},
       W::Vector{Float64}=Vector{Float64}()) = trcorr(M, a, wsg, empt, W)
trcorr(M, a::Vector{uwreal}, W::Vector{Float64},
       wpm::Dict{Int64,Vector{Float64}}) = trcorr(M, a, wsg, wpm, W)
trcorr(M, a::Vector{uwreal}, W::Vector{Float64},
       wpm::Dict{String,Vector{Float64}}) = trcorr(M, a, wsg, dict_name_to_id(wpm), W)
trcorr(M, a::Vector{uwreal},
       wpm::Dict{Int64,Vector{Float64}}) = trcorr(M, a, wsg, wpm, Vector{Float64}())
trcorr(M, a::Vector{uwreal},
       wpm::Dict{String,Vector{Float64}}) = trcorr(M, a, wsg, dict_name_to_id(wpm), Vector{Float64}())


"""
    neid(a::uwreal)

Returns the number of different ensemble ID's contributing to `a`
```@example
using ADerrors # hide
a = uwreal([1.2, 0.2], 12)   # a = 1.2 +/- 0.2
b = uwreal([7.2, 0.5], 13)   # a = 7.2 +/- 0.5

c = a*b
println("Number od ID contributing to c: ", neid(a))
```
"""
neid(a::uwreal)  = ADerrors.unique_ids!(a::uwreal, wsg)

