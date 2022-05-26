include("../src/collections.jl")



col39 = [[1, 2, 4], [4, 5, 7], [1, 7, 8], [1, 2, 5], [4, 5, 8], [2, 7, 8], [1, 2, 8], [2, 4, 5], [5, 7, 8], [2, 5, 8], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7], [6, 7, 8], [7, 8, 9], [1, 8, 9], [1, 2, 9]]
col =  LabelCollection(3,9,col39)

println(blackCliques(col))

wc = whiteCliques(col)
bc = blackCliques(col)

choice = [1,2,5]

white_cliques_w = filter(x -> choice in x, wc)
black_cliques_w = filter(x -> choice in x, bc)

ac_s = union(map(x -> setdiff(choice,intersect(x...)), white_cliques_w)...)
bd_s = union(map(x -> setdiff(union(x...), choice), black_cliques_w)...)


#output: list of crossing pairs [ac,bd], with ac in set1 and bd in set2.
function find_crossing_pairs(set1, set2)
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

# pair : [[a, c], [b, d]]
function pair_is_candidate_for_mutation(collection, mutation_base, pair)
    tmp = setdiff(mutation_base, pair[1])
    if !(union(tmp, [pair[1][1], pair[1][1]]) in collection) ||
    !(union(tmp, [pair[1][1], pair[1][2]]) in collection) ||
    !(union(tmp, [pair[1][2], pair[1][1]]) in collection) ||
    !(union(tmp, [pair[1][2], pair[1][2]]) in collection) ||
        return false
    end
    return true
end

#EXAMPLE
pair_is_candidate_for_mutation(col39, choice, [[5, 1], [4, 8]]) # == true
pair_is_candidate_for_mutation(col39, choice, find_crossing_pairs(ac_s,bd_s)[1]) # == true








