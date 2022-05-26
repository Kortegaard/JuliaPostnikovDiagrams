include("../src/collections.jl")
include("../src/postnikovQuiver.jl")

col39 = [[1, 2, 4], [4, 5, 7], [1, 7, 8], [1, 2, 5], [4, 5, 8], [2, 7, 8], [1, 2, 8], [2, 4, 5], [5, 7, 8], [2, 5, 8], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7], [6, 7, 8], [7, 8, 9], [1, 8, 9], [1, 2, 9]]
col =  LabelCollection(3,9,col39)

function mutations(col::LabelCollection)
    muts = []
    for c in col.collection
        for mut in mutationsAtLabel(col, c)
            push!(muts, mut)
        end
    end
    return muts
end

function mutationsAtLabel(col::LabelCollection, label)
    pairs = findCrossingPairsAtLabel(col, label)
    muts = []
    for pair in pairs
        if !isPairCandidateForMutation(col.collection, label, pair)
            continue 
        end
        
        nCol = union(setdiff(col.collection, [label]), [union(setdiff(label, pair[1]), pair[2])])
        nColObj = LabelCollection(col.k, col.n, nCol)
        sort!(nColObj)
        
        push!(muts, nColObj)
    end
    return muts
end

function findCrossingPairsAtLabel(col::LabelCollection, label)
    wc = whiteCliques(col)
    bc = blackCliques(col)
    
    white_cliques_w = filter(x -> label in x, wc)
    black_cliques_w = filter(x -> label in x, bc)

    ac_s = union(map(x -> setdiff(label,intersect(x...)), white_cliques_w)...)
    bd_s = union(map(x -> setdiff(union(x...), label), black_cliques_w)...)

    return findCrossingPairs(ac_s, bd_s)
end

#output: list of crossing pairs [ac,bd], with ac in set1 and bd in set2.
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

# pair : [[a, c], [b, d]]
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

c36 = [
    LabelCollection(3,6,[[1, 2, 4], [1, 4, 5], [1, 2, 5], [2, 4, 5], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [1, 5, 6], [1, 2, 6]]),
    LabelCollection(3,6,[[1, 2, 4], [1, 4, 5], [1, 3, 4], [1, 4, 6], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [1, 5, 6], [1, 2, 6]])
]

muts = []

for i in c36
    for m in mutations(i)
        if !(m in muts)
            push!(muts, m)
        end
    end
end

for i in 1:100
    mutss = muts
    for i in mutss
        for m in mutations(i)
            if !(m in muts)
                push!(muts, m)
            end
        end
    end
end

for i in 1:length(muts)
    println(isSymmetricCollection(muts[i]))
    #drawPostnikovDiagram(muts[i]; filename="./test/o"*string(i)*".tex", showPostnikovQuiver=true)
    #mycmd = "./test/o"*string(i)*".tex"
    #run(`pdflatex $mycmd -output-directory=./test`)
end
