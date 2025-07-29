using ADerrors, LinearAlgebra

function gen_series!(dt::T where T<:AbstractVector{Float64},τ::Vector{Float64},
                     λ::Vector{Float64};N::Int64=1000,sigma::Float64=1.0)
  nu = zeros(N)
  ap = [t == 0.0 ? 0.0 : exp(-1.0/t) for t in τ]

  for k in eachindex(τ)
    η = randn(N)*sigma
    nu[1] = η[1]
    nu[2:end] = [sqrt(1-ap[k]^2)*η[i] + ap[k]*nu[i-1] for i =2:N]
    [dt[i]+= λ[k]*nu[i] for i in 1:N];
  end
end

@doc raw"""
    gen_cov_series(μ:Vector{Float64}, cov::Matrix{Float64},τ::Vector{Float64},λ::Vector{Float64};N::Int64=1000)

It generate a matrix `length(μ) × N` in which each row is the autocorrelated Montecarlo History of the corresponding μ value.

To generate the MonteCarlo Chain, it first rotate the μ vector into `Y = U^{-1}*μ`, where `U`it the upper triangular matrix obtained by
a Cholesky decomposition, then a Montecarlo chain is generated for each `Y`. Finally, the Montecarlo chain is rotated back with U.
"""
function gen_cov_series(μ::Vector{Float64}, cov::Matrix{Float64},
                        τ::Vector{Float64}, λ::Vector{Float64}; N::Int64=1000)
  np = length(μ)
  res = zeros(np,N); #mc history of each point
  U = LinearAlgebra.cholesky(cov).U;
  Uinv = LinearAlgebra.pinv(U);
  y = Uinv*μ

  res = zeros(np, N)
  for i in 1:np
      gen_series!(view(res,i,1:N),τ,λ,sigma=1.0,N=N)
      res[i,:].+= y[i]
  end

  for t in 1:N
      res[:,t] = U*res[:,t];
  end

  return res;
end


p = let
    mean = [0.1823, 0.2747, 0.3780]
    C = [ 6.4e-7 2.16e-7 4.16e-8;
          2.16e-7 8.1e-7 -3.51e-7;
          4.16e-8 -3.51e-7 1.69e-6]
    series = gen_cov_series(mean,C,[1.0,2.0],[0.3,0.5], N=1_000_000)
    [uwreal(series[i,:], "Epi beta 3.40") for i in 1:3]
end

C = cov(p)
C1 = ADerrors.cov_sym(p)
R = abs.((C-C1)./C).*100

println("Diffence between the two Cov wrt to ADerrors.cov [%]:")
println("\t ",R[1,1], " ",R[1,2], " ",R[1,3])
println("\t ",R[2,1], " ",R[2,2], " ",R[2,3])
println("\t ",R[3,1], " ",R[3,2], " ",R[3,3])

(0==0)
