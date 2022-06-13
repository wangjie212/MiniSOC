function bruteforce(s, n; conf=nothing)
    m = length(s)
    if conf == nothing
        conf = config(n, m)
    end
    for i = 1:length(conf)
        A = zeros(Int16, m+n, n)
        for j = 1:n
            A[conf[i][j][1], j] = 1
            A[conf[i][j][2], j] = 1
            A[j+m, j] = -2
        end
        b = Int16[s;-sum(s);zeros(n-1)]
        if rank(A) == rank([A b])
            return conf[i]
        end
    end
    return nothing
end

function config(n, m)
    conf = Vector{Vector{UInt8}}[]
    t = UInt8[]
    if n > 2
        push!(conf, [[m+2;m+3]])
        push!(t, m+3)
    end
    for i = 1:m
        push!(conf, [[i;m+2]])
        push!(t, m+2)
    end
    for i = 2:n
         nconf = Vector{Vector{UInt8}}[]
         nt = UInt8[]
         for j = 1:length(conf)
             sconf,st = ex_fig(conf[j], i, n, m, t[j])
             append!(nconf, sconf)
             append!(nt, st)
         end
         conf = nconf
         t = nt
     end
     seq1 = [j for j=1:n+m]
     seq2 = [[j for j=1:m]; [j for j=m+2:n+m]]
     red = Int[]
     for i = 1:length(conf)
         if check_nodul(conf[i])
             ind = UInt8[]
             for j = 1:n
                 ind = [ind;conf[i][j]]
             end
             ind = sort(unique(ind))
             if ind == seq1 || ind == seq2
                 push!(red, i)
             end
         end
     end
     return conf[red]
end

function check_nodul(conf)
    for i = 1:length(conf)-1, j = i+1:length(conf)
        if conf[i] == conf[j]
            # || sort([conf[i];m+i]) == sort([conf[j];m+j])
            return false
        end
    end
    return true
end

function ex_fig(conf, i, n, m, et)
    nconf = Vector{Vector{UInt8}}[]
    nt = UInt8[]
    if et >= m + min(i+1, n)
        for j = 1:et-1, k = j+1:et
            if j != m+i && k != m+i
                push!(nconf, [conf; [[j;k]]])
                push!(nt, et)
            end
        end
    end
    if et + 1 <= m + n
        for j = 1:et
            if j != m+i
                push!(nconf, [conf; [[j;et+1]]])
                push!(nt, et+1)
            end
        end
    end
    if et + 2 <= m + n
        push!(nconf, [conf; [[et+1;et+2]]])
        push!(nt, et+2)
    end
    return nconf,nt
end
