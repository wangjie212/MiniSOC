# include("D:\\project\\MiniSOC\\src\\MiniSOC.jl")
using MiniSOC

s = [4;3;2]
conf = GreedyPowertwo(s)

conf = heuristic(s, strategy="GCO")

seq = traversal(s, strategy="GPT")

seq = traversal(s, strategy="GCO")

n = 5
bruteforce(s, n)
