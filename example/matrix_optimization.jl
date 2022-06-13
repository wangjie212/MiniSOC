using JuMP
using MosekTools
using LinearAlgebra
using Dualization

function maxtrace(A, n, s; alg="GPT", QUIET=false)
    # model = Model(optimizer_with_attributes(Mosek.Optimizer))
    model = Model(dual_optimizer(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    w = @variable(model, [1:3], lower_bound=0)
    @constraint(model, sum(w)==1)
    if alg == "GPT"
        conf = GreedyPowertwo(s[1:2])
    else
        conf = fawzi_saunderson(s[2], s[1]+s[2])
    end
    X = Vector{Matrix{VariableRef}}(undef, 2+length(conf))
    for i = 1:2+length(conf)
        X[i] = @variable(model, [1:n^2, 1:n^2])
    end
    @constraint(model, X[1].==kron(w[1]*A[1]+A[2], diagm(ones(n))))
    @constraint(model, X[2].==kron(diagm(ones(n)), w[2]*A[3]+A[4]))
    for i = 1:length(conf)
        M = @variable(model, [1:2n^2, 1:2n^2], PSD)
        @constraint(model, M[1:n^2,1:n^2].==X[conf[i][1]])
        @constraint(model, M[n^2+1:end,n^2+1:end].==X[conf[i][2]])
        @constraint(model, M[1:n^2,n^2+1:end].==X[conf[i][3]])
    end
    if alg == "GPT"
        conf = GreedyPowertwo([s[1]+s[2];s[3]])
    else
        conf = fawzi_saunderson(s[3], s[1]+s[2]+s[3])
    end
    conf = GreedyPowertwo([s[1]+s[2];s[3]])
    Y = Vector{Matrix{VariableRef}}(undef, 2+length(conf))
    for i = 1:2+length(conf)
        Y[i] = @variable(model, [1:n^3, 1:n^3])
    end
    @constraint(model, Y[1].==kron(X[3], diagm(ones(n))))
    @constraint(model, Y[2].==kron(diagm(ones(n^2)), w[3]*A[5]+A[6]))
    for i = 1:length(conf)
        M = @variable(model, [1:2n^3, 1:2n^3], PSD)
        @constraint(model, M[1:n^3,1:n^3].==Y[conf[i][1]])
        @constraint(model, M[n^3+1:end,n^3+1:end].==Y[conf[i][2]])
        @constraint(model, M[1:n^3,n^3+1:end].==Y[conf[i][3]])
    end
    @objective(model, Max, tr(Y[3]))
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if status != MOI.OPTIMAL
       println("termination status: $status")
    end
    println("optimum = $objv")
    return objv
end

function tsallis_entropy(A, n, λ; alg="GPT", QUIET=false)
    # model = Model(optimizer_with_attributes(Mosek.Optimizer))
    model = Model(dual_optimizer(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    w = @variable(model, [1:length(A)], lower_bound=0)
    @constraint(model, sum(w)==1)
    if alg == "GPT"
        conf = GreedyPowertwo([denominator(λ)-numerator(λ);numerator(λ)])
    else
        conf = fawzi_saunderson(numerator(λ), denominator(λ))
    end
    X = Vector{Matrix{Union{VariableRef,Float64}}}(undef, 2+length(conf))
    X[1] = @variable(model, [1:n, 1:n])
    @constraint(model, X[1].==sum(w.*A))
    X[2] = diagm(ones(n))
    for i = 3:2+length(conf)
        X[i] = @variable(model, [1:n, 1:n])
    end
    for i = 1:length(conf)
        M = @variable(model, [1:2n, 1:2n], PSD)
        @constraint(model, M[1:n,1:n].==X[conf[i][1]])
        @constraint(model, M[n+1:end,n+1:end].==X[conf[i][2]])
        @constraint(model, M[1:n,n+1:end].==X[conf[i][3]])
    end
    @variable(model, τ)
    @constraint(model, tr(X[3]-X[1])>=λ*τ)
    @objective(model, Max, τ)
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if status != MOI.OPTIMAL
       println("termination status: $status")
    end
    println("optimum = $objv")
    return objv
end

function fawzi_saunderson(s, t)
    m = UInt8(floor(log(2, t)))
    if t == 2^m
        return tpower(t-s, s, m, 1, 2, 3)
    end
    if 2*s > t
        conf = tpower(t-s, 2^m-(t-s), m, 1, 2, 4)
        k = length(conf) + 4
        append!(conf, tpower(t-2^m, 2^(m+1)-t, m, 2, 3, k))
    else
        conf = tpower(2^m-s, s, m, 1, 2, 4)
        k = length(conf) + 4
        append!(conf, tpower(t-2^m, 2^(m+1)-t, m, 1, 3, k))
    end
    push!(conf, [4;k;3])
    return conf
end

function tpower(s1, s2, m, q1, q2, k)
    g = gcd(s1, 2^m)
    if g > 1
        d = UInt8(log(2, g))
        s1 = Int(s1/2^d)
        s2 = Int(s2/2^d)
        m -= d
    end
    conf = Vector{UInt16}[]
    for i = 1:m-1
        push!(conf, [q1;q2;i+k])
        if s1 < s2
            s2 -= s1
            s1 *= 2
            q1 = i + k
        else
            s1 -= s2
            s2 *= 2
            q2 = i + k
        end
    end
    push!(conf, [q1;q2;k])
    return conf
end

function relative_entropy(ρ, na, nb, λ; QUIET=true)
    n = na * nb
    # model = Model(optimizer_with_attributes(Mosek.Optimizer))
    model = Model(dual_optimizer(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    B = @variable(model, [1:n, 1:n], PSD)
    @constraint(model, tr(B).==1)
    C = @variable(model, [1:n, 1:n], PSD)
    @constraint(model, Tx(B, na, nb).==C)
    conf = GreedyPowertwo([denominator(λ)-numerator(λ);numerator(λ)])
    X = Vector{Matrix{Union{VariableRef,Float64}}}(undef, 2+length(conf))
    X[1] = kron(ρ, diagm(ones(n)))
    for i = 2:2+length(conf)
        X[i] = @variable(model, [1:n^2, 1:n^2])
    end
    @constraint(model, X[2].==kron(diagm(ones(n)), B))
    for i = 1:length(conf)
        M = @variable(model, [1:2n^2, 1:2n^2], PSD)
        @constraint(model, M[1:n^2,1:n^2].==X[conf[i][1]])
        @constraint(model, M[n^2+1:end,n^2+1:end].==X[conf[i][2]])
        @constraint(model, M[1:n^2,n^2+1:end].==X[conf[i][3]])
    end
    @variable(model, τ)
    @variable(model, s)
    @constraint(model, vec(diagm(ones(n)))'*X[3]*vec(diagm(ones(n)))>=s)
    @constraint(model, tr(ρ)-s<=λ*τ)
    @objective(model, Min, τ)
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if status != MOI.OPTIMAL
       println("termination status: $status")
    end
    println("optimum = $objv")
    return objv
end

function Tx(ρ, na, nb)
    A = reshape(ρ, nb, na, nb, na)
    B = Array{VariableRef}(undef, nb, na, nb, na)
    for i = 1:nb
        B[i,:,:,:] = A[:,:,i,:]
    end
    return reshape(B, na*nb, na*nb)
end
