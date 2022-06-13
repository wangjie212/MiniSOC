# MiniSOC
MiniSOC computes (approximately) optimal second-order cone representations for a weighted geometric mean inequality of form $x_1^{\lambda_1}\cdots x_m^{\lambda_m}\ge x_{m+1}$ with $\lambda_i\in\mathbb{Q}^+, i=1,\ldots,m$, $\sum_{i=1}^m\lambda_i=1$. A second-order cone representation is uniquely determined by a set of integer triples $\lbrace(i_k,j_k,t_k)\rbrace_{k=1}^n$ (configuration) so that the weighted geometric mean inequality is equivalent to $x_{i_k}x_{j_k}\ge x_{t_k}^2, k=1,\ldots,n$.  

To use MiniSOC in Julia, run
```Julia
pkg> add https://github.com/wangjie212/MiniSOC
 ```

## Usage
Assume $\lambda_i=\frac{s_i}{p}$, $s_i\in\mathbb{N}$ for $i=1,\ldots,m$ and $p\in\mathbb{N}$. For the functions below, the input $s$ is $(s_1,\ldots,s_m)$.  

Use the fast heuristic algorithm to produce a configuration that determines a
second-order cone representation for $x_1^{\lambda_1}\cdots x_m^{\lambda_m}\ge x_{m+1}$:
```Julia
s = [4;3;2]
conf = GreedyPowertwo(s)
```

Use the brute force algorithm to check whether there exists a second-order cone representation of size $n$:
```Julia
s = [4;3;2]
n = 5
bruteforce(s, n)
```

## References
[Weighted Geometric Mean, Minimum Mediated Set, and Optimal Second Order Cone Representation](), Jie Wang, 2022.  

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn
