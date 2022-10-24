function dectobin(n)
    seq = UInt16[]
    i = 0
    while n > 0
        if isodd(n)
            push!(seq, i)
            n = div(n-1, 2)
        else
            n = div(n, 2)
        end
        i += 1
    end
    return seq
end

function tldeg(a, m)
    if a == 0
        return m
    end
    i = 0
    while iseven(a)
        i += 1
        a = div(a, 2)
    end
    return UInt8(i)
end

function commone(a, b)
    return length(intersect(dectobin(a), dectobin(b)))
end

function heuristic(s; v=length(s)+1, strategy="GPT")
    s = copy(s)
    if v == length(s) + 1
        push!(s, 0)
    end
    p = sum(s)
    m = UInt8(ceil(log(2, p)))
    c = 2^(m-1)
    for i = 1:length(s)-1, j = i+1:length(s)
        if s[i] != 0 && s[i] == s[j]
            if s[i] != c
                push!(s, 2*s[i])
                s[i] = 0
                s[j] = 0
                soc = heuristic(s, v=v, strategy=strategy)
                pushfirst!(soc, [i;j;length(s)])
                return soc
            else
                return [UInt16[i;j;v]]
            end
        end
    end
    r = 2*c - p
    t = argmax(s)
    if 2*s[t] >= p
        if 2*s[t] == p
            s[t] = 0
        elseif s[t] <= c
            s[v] = 2*s[t] - p
            s[t] = 0
        elseif r > 0
            s[v] = r
            s[t] -= c
        else
            s[t] -= c
        end
        soc = heuristic(s, v=length(s)+1, strategy=strategy)
        pushfirst!(soc, [t;length(s)+1;v])
        return soc
    end
    if r > 0
        ds = tldeg.(s, m)
        k = argmin(ds)
        if s[k] <= r && count(i->ds[i]==ds[k], 1:length(s)) == 1
            push!(s, 2*s[k])
            s[k] = 0
            soc = heuristic(s, v=v, strategy=strategy)
            pushfirst!(soc, [k;v;length(s)])
            return soc
        end
        s[v] = r
        k = findfirst(i->s[i]==r&&i!=v, 1:length(s))
        if k != nothing
            push!(s, 2*r)
            s[k] = 0
            s[v] = 0
            soc = heuristic(s, v=v, strategy=strategy)
            pushfirst!(soc, [k;v;length(s)])
            return soc
        end
        p += r
    end
    soc = Vector{UInt16}[]
    for k = 1:1e4
        co = UInt16[]
        cc = Vector{UInt16}[]
        if strategy == "GPT"
            ms = tldeg.(s, m)
            md = minimum(ms)
            for i = 1:length(s)-1, j = i+1:length(s)
                if ms[i] == md && ms[j] == md
                    push!(co, tldeg(s[i]-s[j], m))
                    push!(cc, [i;j])
                end
            end
            mcc = cc[argmax(co)]
            α = min(s[mcc[1]], s[mcc[2]])
        else
            for i = 1:length(s)-1, j = i+1:length(s)
                if s[i] > 0 && s[j] > 0
                    push!(co, commone(s[i], s[j]))
                    push!(cc, [i;j])
                end
            end
            mcc = cc[argmax(co)]
            cmo = intersect(dectobin(s[mcc[1]]), dectobin(s[mcc[2]]))
            α = sum(2 .^cmo)
        end
        push!(s, 2*α)
        s[mcc[1]] -= α
        s[mcc[2]] -= α
        push!(soc, [mcc;length(s)])
        if 2*α >= c
            if 2*α == c
                s[end] = 0
            else
                s[end] -= c
            end
            push!(soc, [length(s);length(s)+1;v])
            append!(soc, heuristic(s, v=length(s)+1, strategy=strategy))
            return soc
        end
        for i = 1:length(s)-1, j = i+1:length(s)
            if s[i] != 0 && s[i] == s[j]
                if s[i] != c
                    push!(s, 2*s[i])
                    s[i] = 0
                    s[j] = 0
                    push!(soc, [i;j;length(s)])
                    append!(soc, heuristic(s, v=v, strategy=strategy))
                else
                    push!(soc, [i;j;v])
                end
                return soc
            end
        end
    end
end

function GreedyPowertwo(s)
    s = copy(s)
    v = length(s) + 1
    w = v
    sub = UInt8[i for i=1:length(s)]
    p = sum(s)
    m = UInt8(ceil(log(2, p)))
    c = 2^(m-1)
    soc = Vector{UInt16}[]
    for iter = 1:1e4
        w += 1
        flag = 0
        for i = 1:length(s)-1, j = i+1:length(s)
            if s[i] == s[j]
                if s[i] != c
                    push!(soc, [sub[i];sub[j];w])
                    s[i] *= 2
                    sub[i] = w
                    deleteat!(s, j)
                    deleteat!(sub, j)
                    flag = 1
                    break
                else
                    push!(soc, [sub[i];sub[j];v])
                    return soc
                end
            end
        end
        if flag == 0
            t = argmax(s)
            r = 2*c - p
            if 2*s[t] >= p
                push!(soc, [sub[t];w;v])
                if 2*s[t] == p
                    p = s[t]
                    deleteat!(s, t)
                    deleteat!(sub, t)
                elseif s[t] <= c
                    temp = s[t]
                    s[t] = 2*s[t] - p
                    sub[t] = v
                    p = temp
                elseif r > 0
                    p = c
                    push!(s, r)
                    push!(sub, v)
                    s[t] -= c
                else
                    p = c
                    s[t] -= c
                end
                v = w
                c = div(c, 2)
                continue
            end
            if r > 0
                ds = tldeg.(s, m)
                k = argmin(ds)
                if s[k] <= r && count(i->ds[i]==ds[k], 1:length(s)) == 1
                    p += s[k]
                    push!(soc, [sub[k];v;w])
                    s[k] *= 2
                    sub[k] = w
                    continue
                end
                k = findfirst(i->s[i]==r, 1:length(s))
                if k != nothing
                    p += r
                    push!(soc, [sub[k];v;w])
                    s[k] = 2*r
                    sub[k] = w
                    continue
                end
                push!(s, r)
                push!(sub, v)
                p += r
            end
            co = UInt8[]
            cc = Vector{UInt8}[]
            ms = tldeg.(s, m)
            md = minimum(ms)
            for i = 1:length(s)-1, j = i+1:length(s)
                if ms[i] == md && ms[j] == md
                    push!(co, tldeg(s[i]-s[j], m))
                    push!(cc, [i;j])
                end
            end
            mcc = cc[argmax(co)]
            α = min(s[mcc[1]], s[mcc[2]])
            if s[mcc[1]] == α
                i1,i2 = mcc[1],mcc[2]
            else
                i1,i2 = mcc[2],mcc[1]
            end
            push!(soc, [sub[i1];sub[i2];w])
            s[i1] *= 2
            s[i2] -= α
            sub[i1] = w
            if 2*α >= c
                w += 1
                if 2*α == c
                    deleteat!(s, i1)
                    deleteat!(sub, i1)
                else
                    s[i1] -= c
                end
                push!(soc, [w-1;w;v])
                v = w
                p = c
                c = div(c, 2)
            end
        end
    end
end

function check_SOC(s, conf)
    m = length(s)
    n = length(conf)
    sort!(conf, lt=(x,y)->isless(x[3],y[3]))
    A = zeros(Int16, m+n, n)
    for j = 1:n
        A[conf[j][1], j] = 1
        A[conf[j][2], j] = 1
        A[j+m, j] = -2
    end
    b = Int16[s;-sum(s);zeros(n-1)]
    if rank(A) == rank([A b])
        return true
    end
    return false
end
