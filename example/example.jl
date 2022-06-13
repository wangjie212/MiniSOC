# include("D:\\project\\MiniSOC\\src\\MiniSOC.jl")
using MiniSOC
using Combinatorics
# using DelimitedFiles
# using AbstractAlgebra

n = 83
part = integer_partitions(n)
m = 8
spart = part[[length(part[i]) == m && gcd(part[i]) == 1 for i=1:length(part)]]
num = Int[]
time = @elapsed begin
for i = 1:length(spart)
    s = spart[i]
    conf = heuristic(s, strategy="GCO")
    push!(num, length(conf))
end
end
println([sum(num), time])

s = Int.(floor.(10^16*rand(8)))
gcd(s) == 1
time = @elapsed begin
conf = GreedyPowertwo(s)
end
println([length(conf), time])

time = @elapsed begin
conf = heuristic(s, strategy="GCO")
end
println([length(conf), time])

# H = hnf(matrix(ZZ, A))

n,m = 6,5
io = open("D://project//MiniSOC//data//config_($n,$m).txt", "w")
conf = config(n, m)
A = zeros(Int, 2*n, length(conf))
for i = 1:length(conf)
    A[:,i] = vcat(conf[i]...)
end
writedlm(io, A, ' ')
close(io)

time = @elapsed begin
seq = traversal(s, strategy="GPT")
end

function brute(s)
    for i = length(s)-1:2*length(s)
        conf = bruteforce(s, i)
        if conf != nothing
            return conf
        end
    end
end

s = [4;4;3]
time = @elapsed begin
io = open("D://project//MiniSOC//data//config_(6,5).txt", "r")
A = readdlm(io, ' ', UInt8)
close(io)
conf = Vector{Vector{Vector{UInt8}}}(undef, size(A,2))
for i = 1:size(A,2)
    conf[i] = [A[j:j+1,i] for j = 1:2:size(A,1)-1]
end
seq = brute(s)
seq = heuristic(s) # 20 0.0002s
end
