#3 state decay problem

#analytical solution
#first state
function populationA(N0, λ, t)
    return N0 * exp(-λ * t)
end

#second state
function populationB(Na0 , λa, λb, t)
    c = (λa * Na0) / (λb - λa)
    return c * (exp(-λa * t) - exp(-λb * t))
end

#third state (final and steady)
function populationC(Na0 , λa, λb, t)
    return Na0 - populationA(Na0, λa, t) - populationB(Na0 , λa, λb, t)
end

#half life to decay constant
#assumes all time units kept consistent (hours for now)
function Hl2Dc(hl)
    return log(2) / hl 
end

#decay constant to half life
Dc2Hl = Hl2Dc

#numerical solution (forward difference aproximation)
#a is isotope up chain, be is isotope being calculated
function d_dt(; λ, λprev, N, Nprev)
    return (λprev * Nprev) - (λ * N)
end

function numAprox(Na0, Nb0, Nc0, λa, λb, dt, tfinal)
    t = collect(0.0:dt:tfinal)
    N = zeros(length(t), 3)
    N[1, :] = [Na0, Nb0, Nc0]

    #[:,1]=A , [:,2]=B , [:,3] = C
    for i in 2:length(t)
        N[i,1] = N[i-1 , 1] + dt* d_dt(λ= λa , λprev=0, N = N[i-1 , 1] , Nprev = 0)
        N[i,2] = N[i-1 , 2] + dt* d_dt(λ= λb , λprev=λa, N = N[i-1 , 2] , Nprev = N[i-1 , 1])
        N[i,3] = N[i-1 , 3] + dt* d_dt(λ = 0, λprev = λb , N = N[i-1,3] , Nprev=N[i-1,2])
    end
    return t , N
end

function t_max_analytical(λa, λb)
    return log(λb / λa) / (λb - λa)
end