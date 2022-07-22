using Random 
using Plots
using Statistics

function get_neighbor(L::Int64)

    N = L*L
    neighbor = zeros(Int64,N,4)

    for site in 1:N  

        if mod(site,L) != 0
            right = site+1
        else
            right = site-L+1
        end    

        if mod(site-1,L) != 0
            left = site-1
        else
            left = site+L-1
        end

        up = site+L

        if up > L*L
            up = up-N
        end    

        down = site-L

        if down < 1
            down = down+N

        end

        neighbor[site,:]=[up,right,down,left]

        
    end

    return neighbor

end


function single_spin_flip(spin::Array{Int64, 2},neighbor::Array{Int64, 2},N::Int64,acceptance::Dict{Int64, Real})

    sum_dE = 0.0
    sum_dM = 0.0

    for flipcount in 1:N

        i = rand(1:N)
        
        j1 = neighbor[i,1]
        j2 = neighbor[i,2]
        j3 = neighbor[i,3]
        j4 = neighbor[i,4]

        dE = 2*spin[i]*(spin[j1] + spin[j2] + spin[j3] + spin[j4])
        

        if rand() < acceptance[dE]

           spin[i] = -spin[i]
           sum_dM += 2*spin[i]
           sum_dE += dE

        end
    end

    return spin, sum_dM, sum_dE

end


function sweep(spin::Array{Int64, 2},neighbor::Array{Int64, 2},n_sweep::Int64,beta::Float64)

    N = length(spin)

    E = zeros(Float64,n_sweep)
    M = zeros(Int64,n_sweep)

    E_sum = zeros(Float64,n_sweep)
    M_sum = zeros(Int64,n_sweep)
    E_2sum = zeros(Float64,n_sweep)
    M_2sum = zeros(Int64,n_sweep)

    acceptance = Dict(-8=>exp(8*beta), -4=>exp(4*beta), 0=>1, 4=>exp(-4*beta), 8=>exp(-8*beta))
    

    E[1] = sum(-spin.*(circshift(spin,(0,1)) .+ circshift(spin,(0,-1)) .+ circshift(spin,(1,0)) .+ circshift(spin,(-1,0))))/2
    M[1] = sum(spin)

    E_sum[1] = E[1]
    M_sum[1] = M[1]

    E_2sum[1] = E[1]^2
    M_2sum[1] = M[1]^2

    for sweepcount in 2:n_sweep

        (spin,sum_dM,sum_dE) = single_spin_flip(spin,neighbor,N,acceptance)
        
        E[sweepcount] = E[sweepcount-1] + sum_dE
        M[sweepcount] = M[sweepcount-1] + sum_dM

        E_sum[sweepcount] = E_sum[sweepcount-1] + E[sweepcount]
        M_sum[sweepcount] = M_sum[sweepcount-1] + M[sweepcount]

        E_2sum[sweepcount] = E_2sum[sweepcount-1] + E[sweepcount]^2
        M_2sum[sweepcount] = M_2sum[sweepcount-1] + M[sweepcount]^2

    end

    return E,M,E_sum,M_sum,E_2sum,M_2sum

end


function main()

    beta = 0.4
    L = 8
    neighbor = get_neighbor(L)
    n_sweep = Int64(1e6)
    spin = ones(Int64,L,L) # cold start

    @time (E,M,E_sum,M_sum,E_2sum,M_2sum) = sweep(spin,neighbor,n_sweep,beta)

    return nothing
end
