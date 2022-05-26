"""
Given a maximal collection `C` [OPS, thm 1.4] describes how and under which circumstances one can mutate a C at a given label.
To check whether we can mutate, we use black and white cliques, which is also explained in [OPS].
"""

include("../src/collections.jl")
include("../src/postnikovQuiver.jl")


"""
    Giving a maximal collection `col` this function calculated a list `muts` of maximal collections of labels which are given by mutating col once.
"""
function mutations(col::LabelCollection)
    muts = []
    for c in col.collection
        for mut in mutationsAtLabel(col, c)
            push!(muts, mut)
        end
    end
    return muts
end


"""
    Given a label `label` in the Labelcollection `col` it computes mutation at that label
"""
function mutationAtLabel(col::LabelCollection, label)
    pairs = findCrossingPairsAtLabel(col, label)
    muts = []
    # TODO: can only mutate if |pairs| = 1 ?
    for pair in pairs
        # Checking if the pair can be mutated at
        if !isPairCandidateForMutation(col.collection, label, pair)
            continue 
        end
        
        # Replacing the label with new label.
        nCol = union(setdiff(col.collection, [label]), [union(setdiff(label, pair[1]), pair[2])])
        nColObj = LabelCollection(col.k, col.n, nCol)
        sort!(nColObj)
        
        push!(muts, nColObj)
    end
    return muts
end


"""
    Given a collection `col` contained a label `label`,
    this function calculates candidates for pairs which can be swapped to use for mutation,
    giving a new maximal non-crossing collection of labels.

    This is done by using black and white cliques,
    to find 'strands' which are near being part of `label` in the sense of being "one away" in the associated postnikov diagram.
"""
function findCrossingPairsAtLabel(col::LabelCollection, label)
    wc = whiteCliques(col)
    bc = blackCliques(col)
    
    # Only keping cliques with label contained
    white_cliques_w = filter(x -> label in x, wc)
    black_cliques_w = filter(x -> label in x, bc)

    # Finds elemets in label which is not in white cliqes.
    ac_s = union(map(x -> setdiff(label,intersect(x...)), white_cliques_w)...)

    # Finds elemets in black cliques which are not in label.
    bd_s = union(map(x -> setdiff(union(x...), label), black_cliques_w)...)

    return findCrossingPairs(ac_s, bd_s)
end


"""
Given two sets of numbers this function calculates pairs [ac,bd] of crossing strands,
i.e. a < b < c < d this respect to a cyclic ordering of Z/nZ.

# Return
    list of crossing pairs [ac,bd], with ac in set1 and bd in set2.
"""
function findCrossingPairs(set1, set2)
    pairs = []
    for ac in combinations(set1,2)
        for bd in combinations(set2,2)
            if isNonCrossing(ac,bd)
                continue
            end
            push!(pairs, [ac,bd])
        end
    end
    return pairs
end


"""
    Checks whether a pair of crossing labels satisfies the condition of [OPS, thm 1.4] for mutations.
"""
function isPairCandidateForMutation(collection, mutation_base, pair)
    tmp = setdiff(mutation_base, pair[1])
    if (!(sort(union(tmp, [pair[1][1], pair[2][1]])) in collection) ||
        !(sort(union(tmp, [pair[1][1], pair[2][2]])) in collection) ||
        !(sort(union(tmp, [pair[1][2], pair[2][1]])) in collection) ||
        !(sort(union(tmp, [pair[1][2], pair[2][2]])) in collection))
        return false
    end
    return true
end
