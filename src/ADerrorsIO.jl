
err(a::uwreal)              = a.err
value(a::uwreal)            = a.mean
derror(a::uwreal)           = a.derr
taui(a::uwreal, i::Int64)   = a.cfd[i].taui
dtaui(a::uwreal, i::Int64)  = a.cfd[i].taui
window(a::uwreal, i::Int64) = a.cfd[i].iw

