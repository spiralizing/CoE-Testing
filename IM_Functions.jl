#Information-based measures
#
"""
    get_tprobs_prev(f_eij, prev_eij, PT_exj)

    Returns a dictionary of probabilites (normalized frequency) from a dictionary
    f_eij and the dictionary of counts prev_eij (previous pool), using a Laplacian
    smoothing within the transition space PT_exij
"""
function get_tprobs_prev(f_eij, prev_eij, PT_exj)
    f_ei = Dict{Any, Int64}()
    for eij in keys(prev_eij)
        ei = split(eij,"+")[1]
        f_ei[ei] = get(f_ei,ei,0) + prev_eij[eij]
    end
    p_ij = Dict{Any,Float64}()
    for eij in keys(f_eij)
        e1 = split(eij,"+")[1]
        if haskey(prev_eij, eij)
            prob_ij = (f_eij[eij] + 1) / (f_ei[e1] + PT_exj[e1])
        else
            prob_ij = 1/PT_exj[e1]
        end
        p_ij[eij] = get(p_ij,eij,prob_ij)
    end
    return p_ij
end
###
"""
    get_tprobs(f_eij)

    Returns probabilities (normalized frequencies) from a dictionary of transition
    counts, the transitions are normalized by the transitions of the first state,
    constructing an stochastic matrix without a matrix.
"""
function get_tprobs(f_eij)
    f_ei = Dict{Any,Int64}()
    for eij in keys(f_eij)
        ei = split(eij,"+")[1]
        f_ei[ei] = get(f_ei,ei,0) + f_eij[eij]
    end

    p_ij = Dict{Any,Float64}()
    for eij in keys(f_eij)
        e1 = split(eij,"+")[1]
        prob_ij = f_eij[eij] / f_ei[e1]
        p_ij[eij] = get(p_ij, eij, prob_ij)
    end
    return p_ij
end
###--
"""
    KLD_Piece_Pool(piece, prev_pool, PT_exj)

    Returns the Kullback-Leibler divergence between a piece (counts of transitions)
    and a pool of previous pieces (counts of transitions).
"""
function KLD_Piece_Pool(piece, prev_pool, PT_exj)
    P_piece = sort(get_tprobs(piece))
    Q_prev = sort(get_tprobs_prev(piece, prev_pool, PT_exj))
    n_e = sum(collect(values(piece)))
    #KL divergence, K(p|q) = Σ p(x)log(p(x)/q(x)P
    P = collect(values(P_piece)) #distribution of the transitions in the piece
    Q = collect(values(Q_prev)) #probabilites of the transitions for the piece in pool
    T = collect(values(piece)) #times the transition is made.
    KLD = mapreduce((x,y,z) -> z*(x* log10(x/y)), +, P, Q, T)
    return KLD
end

###--
"""
    IC_Piece_Pool(piece, prev_pool, PT_exj)

    Returns the information content of a piece (dictionary of pair transitions)
    with the probabilities computed by the MLE in the pool of previous pieces,
    with a laplacian smoothing over the whole transition space.
"""
function IC_Piece_Pool(piece, prev_pool, PT_exj)
    Q_prev = sort(get_tprobs_prev(piece, prev_pool, PT_exj))
    T = collect(values(piece)) #times the transition is made.
    Q = collect(values(Q_prev))
    #information content -Σp(x)
    IC = mapreduce((x,y)-> - y * log10(x), +,Q,T)
    return IC
end

###--
"""
    get_tspace(Xi)

    Returns the number of different transitions x->j in the probability space
    for each x in the alphabet (different states).
"""
function get_tspace(Xi)
    PT_exj = Dict{Any,Int64}() #number of possible transitions in the whole space.
    a_ei = map(x-> split(x,"+")[1],collect(keys(Xi)))
    for i = 1:length(a_ei)
        PT_exj[a_ei[i]] = get(PT_exj, a_ei[i],0)+1
    end
    return PT_exj
end
###############################CODEWORDS
#NOVELTY FUNCTIONS


function get_novelty_measure(piece, Om, N_cws)
    ν_ξ = 0 #the value for the novelty of piece ξ
    ξ = piece
    Ω = collect(keys(Om)) #conventional pool
    Α = collect(setdiff(keys(ξ), keys(Om))) #Novelty pool, elements in the piece that are not in previous works.
    γ = permutedims(hcat(map(x->split(x, "+"), Ω)...)) #re-ordering the pools into 2-dim array.
    α = permutedims(hcat(map(x->split(x, "+"), Α)...))

    s_m = sum(values(ξ)) #size of the piece
    for k in keys(ξ)
        γ_1, γ_2 = split(k,"+") #separate chord 1 and chord 2
        γ_Ω = findall(x->x==γ_1, γ[:,1]) #pairs that start with γ_1 in the Codewordspace
        if isempty(Α)
            s_γ = 0
        else
            γ_a = findall(x->x==γ_1, α[:,1]) #pairs that start with γ_1 in the novelty pool
            s_γ = length(γ_a) #number of pairs that start with γ_1 in novelty pool.
        end
        for e in γ_Ω  #looping over elements that start with γ_1 in the conventional pool
            s_γ += Om[Ω[e]]  #number of pairs in the conventional pool
        end
        s_γ += N_cws
        ix = findall(x->x==k,Ω) #finding the element k in the conventional pool
        if isempty(ix)
            ν_ξ += ξ[k] * log10(s_γ)  #this is the case when the element is a novelty
        else
            ν_ξ += ξ[k] * log10(s_γ / (Om[Ω[ix[1]]] + 1) ) #when the element is in the conventional pool
        end
    end
    return ν_ξ / s_m #dividing by the size of the piece
end

function get_novelty_measure(piece, Om, N_cws; nov_w=1)
    ν_ξ = 0 #the value for the novelty of piece ξ
    ξ = piece
    Ω = collect(keys(Om)) #conventional pool
    Α = collect(setdiff(keys(ξ), keys(Om))) #Novelty pool, elements in the piece that are not in previous works.
    γ = permutedims(hcat(map(x->split(x, "+"), Ω)...)) #re-ordering the pools into 2-dim array.
    α = permutedims(hcat(map(x->split(x, "+"), Α)...))

    s_m = sum(values(ξ)) #size of the piece
    for k in keys(ξ)
        γ_1, γ_2 = split(k,"+") #separate chord 1 and chord 2
        γ_Ω = findall(x->x==γ_1, γ[:,1]) #pairs that start with γ_1 in the Codewordspace
        if isempty(Α)
            s_γ = 0
        else
            γ_a = findall(x->x==γ_1, α[:,1]) #pairs that start with γ_1 in the novelty pool
            s_γ = length(γ_a) #number of pairs that start with γ_1 in novelty pool.
        end
        for e in γ_Ω  #looping over elements that start with γ_1 in the conventional pool
            s_γ += Om[Ω[e]] #number of pairs in the conventional pool
        end
        s_γ += N_cws * nov_w #number of different keys the current key can go to.
        ix = findall(x->x==k,Ω) #finding the element k in the conventional pool
        if isempty(ix)
            ν_ξ += ξ[k]  * log10(s_γ / nov_w)  #this is the case when the element is a novelty
        else
            ν_ξ += ξ[k] * log10((s_γ+nov_w) / (Om[Ω[ix[1]]] + nov_w) ) #when the element is in the conventional pool
        end
    end
    return ν_ξ / s_m #dividing by the size of the piece
end
