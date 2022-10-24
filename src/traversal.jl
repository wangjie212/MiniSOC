function traversal(s; v=length(s)+1, strategy="GPT")
    s = UInt16.(s)
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
                soc = traversal(s, v=v, strategy=strategy)
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
        soc = traversal(s, v=length(s)+1, strategy=strategy)
        pushfirst!(soc, [t;length(s)+1;v])
        return soc
    end
    if r > 0
        ds = tldeg.(s, m)
        k = argmin(ds)
        if s[k] <= r && count(i->ds[i]==ds[k], 1:length(s)) == 1
            push!(s, 2*s[k])
            s[k] = 0
            soc = traversal(s, v=v, strategy=strategy)
            pushfirst!(soc, [k;v;length(s)])
            return soc
        end
        s[v] = r
    end
    index = [s]
    tsoc = [Vector{UInt16}[]]
    tv = UInt16[v]
    for k = 1:1e4
        extsoc = Vector{Vector{UInt16}}[]
        exindex = Vector{UInt16}[]
        extv = UInt16[]
        for j = 1:length(index)
            stsoc,sindex,stv,status = extend(tsoc[j], index[j], m, tv[j], strategy=strategy)
            if status == false
                append!(extsoc, stsoc)
                append!(exindex, sindex)
                append!(extv, stv)
            else
                return stsoc
            end
        end
        tsoc = extsoc
        index = exindex
        tv = extv
    end
end

function extend(soc, s, m, v; strategy="GPT")
    p = sum(s)
    for i = 1:length(s)-1, j = i+1:length(s)
        if s[i] != 0 && s[i] == s[j]
            if 2*s[i] != p
                push!(s, 2*s[i])
                push!(soc, [i;j;length(s)])
                s[i] = 0
                s[j] = 0
                return [soc],[s],[v],false
            else
                push!(soc, [i;j;v])
                return soc,s,v,true
            end
        end
    end
    t = argmax(s)
    if 2*s[t] >= p
        if 2*s[t] == p
            s[t] = 0
        else
            s[t] -= UInt16(p/2)
        end
        push!(s, 0)
        push!(soc, [t;length(s);v])
        return [soc],[s],[length(s)],false
    end
    tsoc = Vector{Vector{UInt16}}[]
    index = Vector{UInt16}[]
    tv = UInt16[]
    lib = Vector{UInt16}[]
    slib = Vector{UInt16}[]
    if strategy == "GPT"
        ms = tldeg.(s, m)
        md = minimum(ms)
        for i = 1:length(s)-1, j = i+1:length(s)
            if ms[i] == md && ms[j] == md
                item = sort([s[i], s[j]])
                if all(l->item!=lib[l], 1:length(lib))
                    push!(lib, item)
                    push!(slib, [i;j])
                end
            end
        end
    else
        cmo = Vector{UInt16}[]
        for i = 1:length(s)-1, j = i+1:length(s)
            if s[i] > 0 && s[j] > 0
                cm = intersect(dectobin(s[i]), dectobin(s[j]))
                if !isempty(cm)
                    item = sort([s[i], s[j]])
                    if all(l->item!=lib[l], 1:length(lib))
                        push!(lib, item)
                        push!(slib, [i;j])
                        push!(cmo, cm)
                    end
                end
            end
        end
    end
    for l = 1:length(lib)
        if strategy == "GPT"
            α = lib[l][1]
        else
            α = sum(2 .^cmo[l])
        end
        cs = copy(s)
        cs[slib[l][1]] -= α
        cs[slib[l][2]] -= α
        push!(cs, 2α)
        push!(index, cs)
        new = push!(copy(soc), [slib[l];length(cs)])
        push!(tsoc, new)
        push!(tv, v)
    end
    return tsoc,index,tv,false
end
