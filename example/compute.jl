# Maximum trace
n = 3
s = [391;280;211]
A = Vector{Matrix{Float64}}(undef, 6)
for i = 1:6
    r = randn(n, n)
    r = r * r'
    A[i] = r/tr(r)
end
@time begin
opt = maxtrace(A, n, s, solver="Mosek", QUIET=true)
end
@time begin
opt = maxtrace(A, n, s, solver="COPT", QUIET=true)
end

# Maximum entropy problem
n = 10
λ = 1//8
A = Vector{Matrix{Float64}}(undef, 10)
for i = 1:10
    r = randn(n, n)
    r = r * r'
    A[i] = r/tr(r)
end
@time begin
opt = tsallis_entropy(A, n, λ, alg="GPT", QUIET=true)
end
@time begin
opt = tsallis_entropy(A, n, λ, alg="FS", QUIET=true)
end

# Tsallis relative entropy
na,nb = 2,2
ρ = randn(na*nb, na*nb)
ρ = ρ * ρ'
ρ = ρ/tr(ρ)
λ = 1//2^8
@time begin
opt = relative_entropy(ρ, na, nb, λ, QUIET=true)
end
