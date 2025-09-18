using LinearAlgebra, ProgressMeter

function σs(i, j)
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

    return coef, spin_flip
end


function Σs(P, s)
    #=
       inputs:
           P: Pauli string
           s: spin config

       returns: tuple with coefficient and resulting spin configuration
    =#
    N = length(s);

    s_ = copy(s)
    coef = 1;
    for i = 1:N
        new_coef, spinflip = σs(P[i], s[i])

        coef *= new_coef
        s_[i] = spinflip
    end

    return (coef, s_)
end

function config_to_index(s)
    N = length(s)
    s_n = (s .+ 1)/2

    s_to_basis = Int(dot(s_n, (2 .^ (0:(N-1))))) + 1

    return s_to_basis
end


function Psi_s(s, vec)
    N = length(s)
    s_n = (s .+ 1)/2

    s_to_basis = Int(dot(s_n, (2 .^ (0:(N-1))))) + 1

    return vec[s_to_basis]
end



function Mn(n, psi; log_base=2)
    @assert n > 1 "Not implemented for n = 1."
    @assert log_base > 0 "Log base must be positive."

    expectation_sum = 0
    N = Int(log2(length(psi)))

    PsiF(s) = Psi_s(s, psi)

    @showprogress for i = 0:(4^N-1)
        P = digits(i, base = 4, pad = N)
        spin_sum = 0
        for j = 0:(2^N-1)
            s = 2 .* digits(j, base = 2, pad = N) .- 1
            alpha, Ps = Σs(P, s)
            spin_sum += conj(PsiF(Ps)) * alpha * psi[j+1]
        end
        expectation_sum += spin_sum^(2*n)
    end

    # default: log base 2 (bits); set log_base = ℯ for natural log (nats)
    out = (1-n)^(-1) * (log(expectation_sum/2^N) / log(log_base))

    if abs(imag(out)) > 10^(-12)
        @warn "Warning imaginary part is non-zero."
        return out
    else
        return real.(out)
    end
end

#n=2 Renyi entropy for mixed states (subtracting log purity)
function M2_mixed(rho; log_base=2)
    expectation_numerator = 0
    expectation_denominator = 0
    N = Int(log2(size(rho)[1]))

    @showprogress for i = 0:(4^N-1)
        P = digits(i, base = 4, pad = N)
        spin_sum = 0
        for j = 0:(2^N-1)
            s = 2 .* digits(j, base = 2, pad = N) .- 1
            alpha, Ps = Σs(P, s)
            spin_sum += rho[j+1, config_to_index(Ps)] * alpha
        end
        expectation_numerator += spin_sum^4
        expectation_denominator += spin_sum^2
    end

    out = -1 * (log(expectation_numerator/expectation_denominator) / log(log_base))

    if abs(imag(out)) > 10^(-12)
        @warn "Warning imaginary part is non-zero."
        return out
    else
        return real(out)
    end

end
