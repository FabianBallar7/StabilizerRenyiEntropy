using Distributed


#Set how many cores you want to run on.
num_procs = 8

addprocs(num_procs)

@everywhere using LinearAlgebra, ProgressMeter, SharedArrays



@everywhere function σs(i,j)
 #=
    inputs:
        i: index of the Pauli. i = 0,1,2,3 => Id, X, Y, Z
        j: index of the spin. j = +-1
    
    returns: tuple with coefficient and resulting spin flip
 =#
    
    if i==0 #Id
        coef = 1
        spin_flip = j 
        
    elseif i==3 #Z
        coef = (+1)*(j==1) + (-1)*(j==-1)
        spin_flip = j 
        
    elseif i==1 #X
        coef = 1
        spin_flip = -1*j
        
    elseif i==2 #Y
        coef = 1im * j
        spin_flip = -1*j
    end
    
    return coef,spin_flip
    
end


@everywhere function Σs(P,s)
 #=
    inputs:
        P: Pauli string
        s: spin config
    
    returns: tuple with coefficient and resulting spin configuration
 =#
    N = length(s);
    
    s_ = copy(s) 
    coef = 1;
    for i=1:N
        new_coef, spinflip = σs(P[i],s[i])
        
        coef *= new_coef
        s_[i] = spinflip
    end
    
    return (coef,s_)
end


@everywhere function config_to_index(s)
    N = length(s)
    s_n = (s .+ 1)/2
    
    s_to_basis = Int(dot(s_n,(2 .^ (0:N-1)))) + 1
    
    return s_to_basis
    
end


@everywhere function Psi_s(s,vec)
    N = length(s)
    s_n = (s .+ 1)/2
    
    s_to_basis = Int(dot(s_n,(2 .^ (0:N-1)))) + 1
    
    return vec[s_to_basis]
    
end



function Mn(n,psi,ncores)
    @assert n > 1 "Not implemented for n = 1."
    
    expectation_sum = 0
    N = Int(log2(length(psi)))
    
    PsiF(s) = Psi_s(s, psi)
    
    out_per_core = SharedArray{ComplexF64}(ncores) 
    
    bin_size = Int(ceil((4^N - 1)/ncores))
    
    worker_batches = collect(Iterators.partition(0:(4^N - 1), bin_size))

    @sync @distributed for k = 1:ncores
        for i in worker_batches[k]
            P = digits(i,base=4,pad=N)
            spin_sum = 0
            for j=0:2^N - 1
                s = 2 .* digits(j,base=2,pad=N) .- 1
                alpha, Ps = Σs(P,s)
                spin_sum += conj(PsiF(Ps)) * alpha * psi[j+1]
            end
            expectation_sum += spin_sum^(2*n)
        end
        out_per_core[k] = expectation_sum
    end    
        
        
    out = (1-n)^(-1) * log(sum(out_per_core)/2^N)
    
        
    if abs(imag(out)) > 10^(-12) 
        @warn "Warning imaginary part is non-zero."
        return out
    else
        return real.(out)
    end
    
end






#n=2 Renyi entropy for mixed states (subtracting log purity)
function M2_mixed(rho,ncores)
    
    expectation_numerator = 0
    expectation_denominator = 0
    N = Int(log2(size(rho)[1]))
    
    out_per_core_num = SharedArray{ComplexF64}(ncores) 
    out_per_core_den = SharedArray{ComplexF64}(ncores) 
    
    
    bin_size = Int(ceil((4^N - 1)/ncores))
    
    worker_batches = collect(Iterators.partition(0:(4^N - 1), bin_size))

    @sync @distributed for k = 1:ncores
        for i in worker_batches[k]
            P = digits(i,base=4,pad=N)
            spin_sum = 0
            for j=0:2^N - 1
                s = 2 .* digits(j,base=2,pad=N) .- 1
                alpha, Ps = Σs(P,s)
                spin_sum += rho[j+1,config_to_index(Ps)] * alpha
            end
            expectation_numerator += spin_sum^4
            expectation_denominator += spin_sum^2
        end
        out_per_core_num[k] = expectation_numerator
        out_per_core_den[k] = expectation_denominator
    end    

    out = -1*log(sum(out_per_core_num)/sum(out_per_core_den))
    
    if abs(imag(out)) > 10^(-12) 
        @warn "Warning imaginary part is non-zero."
        return out
    else
        return real(out)
    end
    
end















