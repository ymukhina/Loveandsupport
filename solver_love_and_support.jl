using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff 
using Random

const Ptype = QQMPolyRingElem

# -------- estimate support for f_min based on Theorems 1 and 2 -------- #

function f_min_support_bivariate(ode::ODE)
    @assert length(ode.x_vars) == 2 "system has more than two variables."
    x1, x2 = sort(ode.x_vars, rev = true)
    g1 = ode.x_equations[x1]
    g2 = ode.x_equations[x2]
    d1 = total_degree(g1)
    d2 = total_degree(g2)
    if d1 <= d2
        A = [1 d1 (d1 + d2 - 1); -1 0 0; 0 -1 0; 0 0 -1]
        b = [d1 * (d1 + d2 - 1), 0, 0, 0]
    else
        A = [1 d1 (2 * d1 - 1); 1 d2 (d1 + d2 - 1); -1 0 0; 0 -1 0; 0 0 -1]
        b = [d1 * (2 * d1 - 1), d1 * (d1 + d2 - 1), 0, 0, 0]
    end
    return collect(lattice_points(polyhedron(A, b)))
end

function f_min_support(ode::ODE, jacobian_rank::Int; info = true)
    # I guess by sorting you try to ensure that the first x is the output, i am not sure
    # this would always work. You could just choose x such that the only y_equation is equal to x
    x = sort(ode.x_vars, rev = true)
    info && @info "Inferred order $x"
    n = jacobian_rank
    gs = [ode.x_equations[xi] for xi in x]
    d1 = total_degree(gs[1])
    @assert d1 > 0 "d1 = 0"
    D = maximum(total_degree, gs[2:end])
    D = max(D, 0)
    info && @info "We have d1 = $d1 and D = $D"
    if d1 <= D
        ineq_lhs = reshape([1, [d1 + (k - 1) * (D - 1) for k in 1:n]...], 1, n + 1)
        ineq_rhs = [prod([d1 + (k - 1) * (D - 1) for k in 1:n])]
        A = vcat(matrix(QQ, ineq_lhs), -identity_matrix(QQ, n + 1))
        b = vcat(ineq_rhs, zeros(QQ, n + 1))
    else
        ineq_lhs1 = [k <= l ? k * (D - 1) + 1 : 0 for l in 0:(n - 1), k in 0:n]
        ineq_lhs2 = zeros(Int, n, n + 1)
        for l in 0:(n - 1)
            for i in 1:(n - l)
                ineq_lhs2[l + 1, i + l + 1] = i * (d1 - 1) + l * (D - 1) + 1
            end
        end
        # Split maybe?
        ineq_rhs = [prod(Vector{Int}([d1 + (k - 1) * (D - 1) for k in 1:l])) * prod(Vector{Int}([i * (d1 - 1) + l * (D - 1) + 1 for i in 1:(n - l)]))
                    for l in 0:(n - 1)]
        A = vcat(matrix(QQ, ineq_lhs1 + ineq_lhs2), -identity_matrix(QQ, n + 1))
        b = vcat(ineq_rhs, zeros(QQ, n + 1))
    end
    return collect(lattice_points(polyhedron(A, b)))
end

# ---------------------------------------------------------------------- #


# -------- Compute f_min using an ansatz equation -------- #
    
function build_matrix_power_series(F, ode, support; info = true)

    n = length(ode.x_vars)
    y = first(ode.y_vars)
    nterms = length(support) + n
        
    # random initial conditions
    ic = Dict([v => rand(F) for v in ode.x_vars]...)
    # no parameters, no inputs
    par = empty(ic)
    inp = empty(Dict(first(ode.x_vars) => [one(F)]))
           
     
      
    ps_soltime = @elapsed ps_sol = power_series_solution(ode, par, ic, inp, nterms)
    info && @info "Power series solution computed in $ps_soltime"

    start_system_time = time()
    pss = [ps_sol[y]]
    for i in 1:n
        push!(pss, ps_diff(pss[end]))
    end
            
    sort!(support, by = sum)
            
    prods = eval_at_support(support, pss)
                
    ls = matrix_eff(prods, F, nterms - n) 
           
    info && @info "System created in $(time() - start_system_time)"
    return ls
end
                    
                    
function check(pol, ode::ODE)
                        
    x = first(sort(ode.x_vars, rev = true))                     
    n = length(ode.x_vars)
    
    dervs = var_derivatives(n, ode, x) 

    res = pol(dervs...) 
                        
   return iszero(res), res
end
                        
function deg_after_substitution(pol, degs)
           
      mons = collect(Oscar.exponents(pol))                               
      k = findmax([sum(mons[j] .* degs) for j in 1:length(mons)])[1]
                                                                  
      return(k)                          
end                                
                        
                        
function ShZ_check(pol, ode::ODE, target_prob) 
                            
    x = first(sort(ode.x_vars, rev = true))                     
    n = length(ode.x_vars)
    jac = jacobian_check(ode) 
                                                                             
    dervs = var_derivatives(jac, ode, x) 
       
    D = [total_degree(d) for d in dervs]                        
                                
    deg_bnd = deg_after_substitution(pol, D)                                                                       
    
    N = Int(ceil(deg_bnd / (1 - target_prob)))                               
                                
    vec = [rand(1:N) for _ in 1:(n + 1)]
     
    evals = [derv(vec...) for derv in dervs]                        
                       
    res = pol(evals...)
                            
   return iszero(res)
end                            

function build_matrix_multipoint(F, ode, support; info = true)
                                
    x = first(sort(ode.x_vars, rev = true))
    n = length(ode.x_vars)
    @info "computing derivatives"
    dervs = var_derivatives(n, ode, x)
    @info "done"               
    
    sort!(support, by = sum)
    
    # reorganize to compute efficiently and in-place
    lsup = length(support) # N
    # TODO: create a Nemo matrix (with matrix_space) right away
    mat = Array{elem_type(F), 2}(undef, lsup, lsup)
        for i in 1:lsup # N iterations
            vec = [rand(F) for _ in 1:(n + 1)]
            evals = [derv(vec...) for derv in dervs]
            #prods = [prod(evals .^ exp) for exp in support]
            #push!(prods, [prod(evals .^ exp) for exp in support])
            #evsup = eval_at_support(support, evals)
                for j in 1:lsup # N iterations
                    mat[i, j] = prod([evals[k]^support[j][k] for k in 1:(n+1)]) # O(d) , where d - degree of the minimal polynomial
                    #mat[i, j] = prod(evals .^ support[j])
                end
        end # O(N^2 * d)

    # print runtime
    return matrix(mat)
end
                                        
# return the order of the largest upper left non vanishing minor of the jacobian matrix                          
function jacobian_check(ode)
                                            
    x = first(sort(ode.x_vars, rev = true))                                          
    n = length(ode.x_vars)  
                    
    dervs = var_derivatives(n-1, ode, x)
    J = jacobian_matrix(dervs)[1:n, :]
    res = rank(J)                                        
 
    return res                                                                                              
end                                            

function build_matrix_multipoint_fast(F, ode, jac, support; info = true)
    x = first(sort(ode.x_vars, rev = true))
    n = length(ode.x_vars)
    @info "computing derivatives"
    dervs = var_derivatives(jac, ode, x)                                      
    @info "done"               

    support = [Vector{Int64}(p) for p in support]
    # sort!(support, by = sum)
    
    lsup = length(support)                                                    
    S = matrix_space(F, lsup, lsup)
    M = zero(S)
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    # filling the columns corresponding to the derivatives
        for i in 1:lsup
            M[i, 1] = 1
            vec = [rand(F) for _ in 1:(n + 1) ]                                      
            evals = [derv(vec...) for derv in dervs]
                                      
                                                        
                for j in 1:(jac + 1)
                    supp = [(k == j) ? 1 : 0 for k in 1: (jac + 1) ]
                    ind = supp_to_index[supp]
                    M[i, ind] = evals[j]
                end
        end

    # filling the rest of the columns
        for i in (jac + 3):lsup
            supp = support[i]
            nonzero_ind = findfirst(x -> x > 0, supp)
            var_ind = supp_to_index[[(k == nonzero_ind) ? 1 : 0 for k in 1: (jac + 1) ]]
            pred_supp = copy(supp)
            pred_supp[nonzero_ind] -= 1
            pred_ind = supp_to_index[pred_supp]
                for j in 1:lsup
                    M[j, i] = M[j, pred_ind] * M[j, var_ind]
                end
        end
      
  return M
end


function eliminate_with_love_and_support_modp(ode::ODE, p::Int; info = true)
                                                                    
    @assert is_probable_prime(p) "This is not a prime number, Yulia!"

    ode_mod_p = reduce_ode_mod_p(ode, p)
    x = sort(ode_mod_p.x_vars, rev = true)
    n = length(x)
    F = Nemo.Native.GF(p)     

    # take support as an argument and compute only once
    # compute Newton polytope of f_min
    jac = jacobian_check(ode) 
                                                                        
    possible_supp = f_min_support(ode, jac)
                                                               
                                                              
    l = length(possible_supp)
    info && @info "The size of the estimates support is $(length(possible_supp))"
    #@info possible_supp

    # make an extra argument to choose a function
    # tim1 = @elapsed ls = build_matrix_power_series(F, ode_mod_p, possible_supp, info = info)
    sort!(possible_supp, by = sum)
    tim2 = @elapsed ls = build_matrix_multipoint_fast(F, ode_mod_p, jac, possible_supp, info = info)
    
                                                                    
    info && @info "eval method $(tim2)"

    info && @info "linear system dims $(size(ls))"
    
    system_soltime = @elapsed ker = kernel(ls, side=:right)
    info && @info "Linear system solved in $system_soltime"
            dim = size(ker)[2]
    info && @info "The dimension of the solution space is $(dim)"

    start_constructing_time = time()

    R, _ = polynomial_ring(F, ["x1", ["x1^($i)" for i in 1:jac]...])
            
    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]    
    
    g = gcd([sum([s * m for (s, m) in zip(ker[:, i], mons)]) for i in 1:dim])
       
    info && @info "The resulting polynomial computes in $(time() - start_constructing_time)"

    return g
end
        

# -------------------------------------------------------- #

# -------- Faster evaluation of many monomials in power series -------- #
        
function matrix_eff(vals, F, N)
                                                                        
    S = matrix_space(F, N, N)
    ls = zero(S)
        for i in 1:N
            for j in 1:N
                ls[j,i] = coeff(vals[i], j-1)
            end
        end
                                                                                
   return ls
end

#friendly helper
function one_thing(e, vals, cacher)
    if !haskey(cacher, e)
        for i in 1:length(e)
            if e[i] != 0 #[1,?]
                e_love = copy(e)
                e_love[i] -= 1
                res = vals[i] * one_thing(e_love, vals, cacher)
                cacher[e] = res
            end
         end  
    end
  return cacher[e]
end    

function eval_at_support(supp, vals)
    s = sort(supp, by = sum) 
    cacher = Dict(s[1] => one(parent(first(vals))))
            for e in s[2:end] #index
                one_thing(e, vals, cacher)                 
            end
return [cacher[s] for s in supp]
end
    
  
function qq_to_mod(a::QQFieldElem, p)
  return numerator(a) * invmod(denominator(a), ZZ(p))
end
                                                                                                                
function final(ode::ODE) 
                                                                                                                
    g, a = eliminate_with_love_and_support(ode, 2^17-1)
                                                                                                                    
         while ShZ_check(g, ode, 0.9) == false
                 println(total_degree(g-1))                                                                                       
                 a = Hecke.next_prime(a)
                 g, a = eliminate_with_love_and_support(ode, a)
                 println(a) 
         end
                                                                                                                         
return g                                                                                                             
end                                                                                                                    

function eliminate_with_love_and_support(ode::ODE, a::Int)
   jac = jacobian_check(ode) 
   possible_supp = f_min_support(ode, jac)
   l_supp = length(possible_supp)
   #x = sort(ode.x_vars, rev = true)
   #n = length(x)
   R, _ = polynomial_ring(QQ, ["x1", ["x1^($i)" for i in 1:jac]...]) 
   sort!(possible_supp, by = exp -> prod(gens(R) .^ exp), rev = true)
   mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]

   sol_vector = zeros(QQ, l_supp)
   crts = zeros(ZZ, l_supp)
   found_cand = falses(l_supp)
   is_stable = falses(l_supp)

   prod_of_done_primes = one(ZZ)
   prim_cnt = 0
                                                                                                                                                                                                                            
   while !all(is_stable)
       @label nxt_prm                                                                                                                 
       a = Hecke.next_prime(a)
       p = ZZ(a)                                        
       prim_cnt += 1
                                                                                                                       
       l_succ = length(findall(is_stable))
       @info "Chose $prim_cnt th prime $p, $(l_succ) stable coefficients"
       sol_mod_p = eliminate_with_love_and_support_modp(ode, Int(p))
       sol_vector_mod_p = [coeff(sol_mod_p, Vector{Int}(exp)) for exp in possible_supp]

       for (i, a) in enumerate(sol_vector_mod_p)
           is_stable[i] && continue

           if found_cand[i]
               if divides(denominator(sol_vector[i]), p)[1]
                   @info "bad prime, restarting"
                   @goto nxt_prm
               end
               sol_i_mod_p = qq_to_mod(sol_vector[i], p)
               if sol_i_mod_p == sol_vector_mod_p[i]
                   is_stable[i] = true
                   continue
               end
           end

           crts[i] = crt(lift(ZZ, a), p, ZZ(crts[i]), prod_of_done_primes)
           
           try
               succ, r, s = rational_reconstruction(crts[i],
                                                    ZZ(p*prod_of_done_primes))
               if succ
                   sol_vector[i] = r//s
                   found_cand[i] = true
               end
                catch
                    @info "flint crisis"
                end
       end

       prod_of_done_primes *= p
   end 
    
    sort!(mons, rev = true)
    g = sum([s * m for (s, m) in zip(sol_vector, mons)])
   return g, a
end 
        
# -------- Helpers for interpolation Ansatz -------- #
        
function Lie_derivative(pol, ode)
    result = zero(pol)
        for v in vars(pol)
            result += derivative(pol, v) * ode.x_equations[v]
        end
   return result
end
      
        
function var_derivatives(ord, ode, var)
    result = [var]
            
        for i in 1:ord
            push!(result, Lie_derivative(last(result), ode))
        end
                
   return result
end
            
function evaluate_exp_vector(vec, vals)
    return prod(vals .^ vec)
end

# ————————————————— #
