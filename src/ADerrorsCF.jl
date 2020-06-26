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
const DO_BIN = false

function get_nbin_vec(nd::Array{Int64, 1})::Int64
    
    nbin::Int64 = 1
    if (!DO_BIN)
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

function uwcls(data::Vector{Float64}, id::Int64, ws::wspace, iv::Vector{Int64})
    if (length(data) == 2)
        ws.nob += 1
        new  = fbd(1, 1, [data[2]], [1], Dict{Int64,Vector{Complex{Float64}}}())
        push!(ws.fluc, new)
        push!(ws.map_nob, id)
        if (!haskey(ws.map_ids, id))
            ws.map_ids[id] = ws.nob
        end
        
        p = [false for n in 1:ws.nob]
        p[end] = true
        d = [0.0 for n in 1:ws.nob]
        d[end] = 1.0
        return uwreal(data[1], 0.0, 0.0,
                      p, d, Vector{Int64}(), Vector{cfdata}())
    else
        if (sum(iv) != length(data))
            ArgumentError("Sum of replica length does not match number of measurements")
        end
        ws.nob += 1
        avg = Statistics.mean(data)
        push!(ws.map_nob, id)
        if (!haskey(ws.map_ids, id))
            ws.map_ids[id] = ws.nob
        end

        fseries = Dict{Int64,Vector{Complex{Float64}}}()
        nbin = get_nbin_vec(iv)
        is = 1
        for i in 1:length(iv)
            ie = is + iv[i] - 1
            nbdt = div(iv[i], nbin)
            datapad = [bin_data(data[is:ie] .- avg, nbin);
                       zeros(Float64, nextpow(2,2*nbdt+1)-nbdt) ]
            fseries[i] = FFTW.fft(datapad)
            is = ie + 1
        end

        new = fbd(length(data), nbin,
                  data .- avg,
                  iv,
                  fseries)
        push!(ws.fluc, new)
        
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

    avg = Statistics.mean(data)
    dt = fill(avg, nms)
    for n in 1:length(idm)
        dt[idm[n]] = data[n] - avg
    end
    dt .= (nms/length(idm)) .* dt
    
    ws.nob += 1
    push!(ws.map_nob, id)
    if (!haskey(ws.map_ids, id))
        ws.map_ids[id] = ws.nob
    end

    fseries = Dict{Int64,Vector{Complex{Float64}}}()
    nbin = get_nbin_vec(iv)
    is = 1
    for i in 1:length(iv)
        ie = is + iv[i] - 1
        nbdt = div(iv[i], nbin)
        datapad = [bin_data(dt[is:ie], nbin);
                   zeros(Float64, nextpow(2,2*nbdt+1)-nbdt) ]
        fseries[i] = FFTW.fft(datapad)
        is = ie + 1
    end

    new = fbd(nms, nbin,
              dt,
              iv,
              fseries)
    push!(ws.fluc, new)

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
                if (any(a.ids[1:end] == ws.map_nob[i]))
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
                return i
            else
                tau = stau/log((2.0*tiw+1.0)/(2.0*tiw-1.0))
                gw = exp(-(i-1.0)/tau) - tau/sqrt((i-1.0)*nd)
                if (gw < 0.0)
                    return i
                end
            end
        end
    end
    
end

function uwerror(a::uwreal, ws::wspace, wpm::Dict{Int64,Vector{Float64}})

    nid = unique_ids!(a, ws)

    if (length(a.cfd) == 0)
        for j in 1:nid
            new = cfdata(0.0, 0.0, 0.0, 0, Vector{Float64}(), Vector{Float64}())
            idx = ws.map_ids[a.ids[j]]
            nd  = ws.fluc[idx].nd
            if (nd != 1)
                nd_eff = div(nd, ws.fluc[idx].ibn)
                (nt,ip) = findmax(ws.fluc[idx].ivrep)
                nt = div(nt, 2*ws.fluc[idx].ibn)
                nrep   = length(ws.fluc[idx].ivrep)
                ftemp  = Dict{Int64,Vector{Complex{Float64}}}()
                for i in 1:nrep
                    ns = length(ws.fluc[idx].fourier[i])
                    ftemp[i] = zeros(Complex{Float64}, ns)
                end
            end
            
            for i in 1:length(a.prop)
                if (a.prop[i] && (ws.map_nob[i] == a.ids[j]))
                    if (nd == 1)
                        new.var = new.var + a.der[i]*ws.fluc[i].delta[1]
                    else
                        for k in 1:nrep
                            ftemp[k] = ftemp[k] + a.der[i]*ws.fluc[i].fourier[k]
                        end
                    end
                end
            end

            if (nd == 1)
                new.var = new.var^2
            else
                new.gamm = zeros(nt)
                for k in 1:nrep
                    ftemp[k] .= ftemp[k].*conj(ftemp[k])
                    FFTW.ifft!(ftemp[k])

                    for ig in 1:min(nt,length(ftemp[k]))
                        new.gamm[ig] = new.gamm[ig] + real(ftemp[k][ig])
                    end
                end
                    for ig in 1:nt
                    nrcnt = count(map(x -> div(x, 2*ws.fluc[idx].ibn), ws.fluc[idx].ivrep) .> ig-1)
                    new.gamm[ig] = new.gamm[ig] / (nd_eff - nrcnt*(ig-1))
                end

                iw = wopt_ulli(nd_eff, DEFAULT_STAU, new.gamm)
                new.iw = iw
                
                dbias = new.gamm[1] + 2.0*sum(new.gamm[2:iw])
                new.gamm .= new.gamm .+ dbias/nd_eff

                new.drho = zeros(nt)
                if (new.gamm[1] != 0.0)
                    for i in 1:nt
                        is = max(1, i-iw-2) + 1
                        ie = i + iw-1
                        for k in is:ie
                            if (k < nt+1)
                                cont = -2.0*new.gamm[k]*new.gamm[i]/new.gamm[1]^2
                            else
                                cont = 0.0
                            end

                            if ((i+k-2) < nt)
                                cont = cont + new.gamm[i+k-1]/new.gamm[1]
                            end
                            if (abs(i-k) < nt)
                                cont = cont + new.gamm[abs(i-k)+1]/new.gamm[1]
                            end
                            new.drho[i] = new.drho[i] + cont^2
                        end
                        new.drho[i] = sqrt(new.drho[i]/nd_eff)
                    end
                else
                    new.var = new.var^2
                end
            end

            push!(a.cfd, new)
        end
    end

    
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
            ibn     = ws.fluc[idx].ibn
            nd_eff  = div(ws.fluc[idx].nd, ibn)
            (nt,ip) = findmax(ws.fluc[idx].ivrep)
            nt = div(nt, 2*ibn)
            if (a.cfd[j].gamm[1] == 0.0)
                a.cfd[j].taui  = 0.0
                a.cfd[j].dtaui = 0.0
                vti = 0.0
                iw  = 0
            else
                wp = zeros(4)
                wp = get(wpm, a.ids[j], [-1.0,-1.0,-1.0,-1.0])
                if (wp[1] > 0.0)
                    a.cfd[j].iw = round(wp[1])
                elseif (wp[2] > 0.0)
                    a.cfd[j].iw = wopt_ulli(nd_eff, wp[2], a.cfd[j].gamm)
                end
                
                if (wp[3] > 0.0)
                    for k in 2:nt
                        if (a.cfd[j].drho[k]*wp[3] > a.cfd[j].gamm[k]/a.cfd[j].gamm[1])
                            continue
                        end
                    end
                    a.cfd[j].iw = k-1
                end

                if (wp[4] > 0.0)
                    texp = wp[4]/ibn
                else
                    texp = 0.0
                end
                iw = a.cfd[j].iw
                
                vti = 0.5 + sum(a.cfd[j].gamm[2:iw])/a.cfd[j].gamm[1]
                a.cfd[j].dtaui = sqrt(vti^2 * (4.0*iw-2.0*vti+2.0)/nd_eff)
                a.cfd[j].taui  = vti + texp*a.cfd[j].gamm[iw+1]/a.cfd[j].gamm[1]
            end
            
            a.cfd[j].var = a.cfd[j].gamm[1] * 2.0*a.cfd[j].taui/nd_eff
        end

        a.err = a.err + a.cfd[j].var
        if (iw > 1)
            a.derr = a.derr + vti*(a.cfd[j].iw-0.5)/nd_eff
        end
    end
    
    a.err  = sqrt(a.err)
    a.derr = sqrt(a.derr)
end

function unique_ids_multi(a::Vector{uwreal})

    is = 0
    for i in 1:length(a)
        is = is + unique_ids!(a)
    end

    n = 0
    ids = Vector{Int64}()
    for k in 1:length(a)
        for i in 1:length(a[k].prop)
            if (a[k].prop[i])
                if (any(a.ids[1:end] == ws.map_nob[i]))
                    continue
                end
                n = n + 1
                push!(ids, map_nob[i])
            end
        end
    end

    return ids
end


function cov(a::Vector{uwreal}, ws::wspace, wpm::Dict{Int64,Vector{Float64}})

    for i in 1:length(a)
        uwerr(a, wpm)
    end
    ids = unique_ids_multi(a)
    nid = length(ids)
    iw  = zeros(Int64, nid)

    for k in 1:length(a)
        for j in 1:length(a[k].ids)
            for i in 1:nid
                if (a[k].ids[j] == ids[j])
                    iw[i] = max(iw[i], a[k].cfd[j].iw)
                end
            end
        end
    end

    wopt = Dict{Int64,Vector{Float64}}()
    for i in 1:nid
        wopt[ids[i]] = [convert(Float64, iw[i]), -1.0, -1.0, -1.0]
    end

    cov = zeros(Float64, length(a), length(a))
    for k in 1:length(a)
        uwerr(a[k], wopt)
        cov[k,k] = a[k].err^2
    end

    for j in 1:nid
        idx  = ws.map_ids[a.ids[j]]
        nd   = ws.fluc[idx].nd
        nrep = length(ws.fluc[idx].ivrep)
        for k1 in 1:length(a)-1
            v1 = 0.0
            for i in 1:length(a[k1].prop)
                if (a[k1].prop[i] && (ws.map_nob[i] == ids[j]))
                    if (nd == 1)
                        v1 = v1 + a[k1].der[i]*ws.fluc[i].delta[1]
                    else
                        # AQUI
                    end
                end
            end
            for k2 in k1+1:length(a)
                v2 = 0.0
                for i in 1:length(a[k2].prop)
                    if (a[k2].prop[i] && (ws.map_nob[i] == ids[j]))
                        if (nd == 1)
                            v2 = v2 + a[k2].der[i]*ws.fluc[i].delta[1]
                        end
                    end
                end
            end

            if (nd == 1)
                cov[k1,k2] = cov[k1,k2] + v1*v2
            end
        end
    end
    
    
end


##                                         ##
# Module workspace and function definitions #
##                                         ##


wsg = ADerrors.wspace(similar(Vector{ADerrors.fbd}, 0),
                0,
                similar(Vector{Int64}, 0),
                Dict{Int64, Int64}())

empt = Dict{Int64,Vector{Float64}}()

uwreal(data::Vector{Float64}, id::Int64) = ADerrors.uwcls(data::Vector{Float64}, id::Int64, wsg, [length(data)])
uwreal(data::Vector{Float64}, id::Int64, iv::Vector{Int64}) = ADerrors.uwcls(data::Vector{Float64}, id::Int64, wsg, iv)

uwerr(a::uwreal) = ADerrors.uwerror(a::uwreal, wsg, empt)
uwerr(a::uwreal, wpm::Dict{Int64,Vector{Float64}}) = ADerrors.uwerror(a::uwreal, wsg, wpm::Dict{Int64,Vector{Float64}})
neid(a::uwreal)  = ADerrors.unique_ids!(a::uwreal, wsg)

