using Combinatorics

mutable struct LabelCollection
    k::Int
    n::Int

    collection::Vector{Vector{Int}}
end

Base.copy(s::LabelCollection) = LabelCollection(s.k, s.n, s.collection)
Base.hash(c::LabelCollection) = Base.hash(c.collection)
Base.:(==)(a::LabelCollection, b::LabelCollection) = (a.k == b.k && a.n == b.n && isEquivalentToCollectionUpToRotation(a.n, a.collection, b.collection))

"""
    So
"""
function sort!(col::LabelCollection)
    col.collection = sort(map( x->sort(x), col.collection))
end

"""
    true if and only if a<b<c cyclically
"""
function isCyclicOrdered(a,b,c)
    diff_b = b-a;
    diff_c = c-a;
    if diff_b > 0 && diff_c > 0
        return diff_b < diff_c
    elseif diff_b < 0 && diff_c < 0
        return diff_b < diff_c
    elseif diff_c < 0
        return true
    end
    return false
end

"""
    Checks if two labels a and b is a non-crossing pair of labels.
"""
function isNonCrossing(a,b)
    aa = setdiff(a,b);
    bb = setdiff(b,a);
    
    a_index = 1
    b_index = 1

    diff_length = 1;
    curr_obj = min(aa[1],bb[1]);
    at_a = argmin([aa[1],bb[1]]) == 1;
    if at_a; a_index = (a_index +1 > length(aa) ? 1 : a_index +1) end
    if !at_a; b_index = (b_index +1 > length(bb) ? 1 : b_index +1) end

    for i in 1:length(aa) + length(bb) + 4
        if isCyclicOrdered(curr_obj, aa[a_index], bb[b_index])
            diff_length = (at_a ? 0 : diff_length+1)
            at_a = true;
            curr_obj = aa[a_index]
            a_index = (a_index +1 > length(aa) ? 1 : a_index +1)
        else
            diff_length = (!at_a ? 0 : diff_length+1)
            at_a = false;
            curr_obj = bb[b_index]
            b_index = (b_index +1 > length(bb) ? 1 : b_index +1)
        end
        if diff_length >= 4; return false end
    end
    return true
end

"""
    Checks if a collection of labels is non-crossing
"""
function isNonCrossingCollection(list)
    for i in 1:length(list)
        for j in i+1:length(list)
            if i == j; continue; end
            if !isNonCrossing(list[i],list[j]); return false; end
        end
    end
    return true
end

function isNonCrossingCollection(col::LabelCollection)
    return isNonCrossingCollection(col.collection)
end


"""
    Checking if a collection of Noncorssing labels are maximal
    see [PTZ ,thm. 1.8], [OPS thm. 4.7]
"""
function isMaximalNonCrossingCollection(k, n, list)
    return isNonCrossingCollection(list) && length(list) == k*(n-k)+1
    #if !isNonCrossingCollection(list) return false end
    #for it in combinations(1:n,k)
    #    if it in list; continue; end
    #    if isNonCrossingCollection(union(list, [it])) return false end
    #end
    #return true;
end


"""
    Checking if a collection of Noncorssing labels are maximal
"""
function isMaximalNonCrossingCollection(k, n, list, exclude_projectives)
    if exclude_projectives && length(list) == k*(n-k)-n+1; return isNonCrossingCollection(list); end;
    if exclude_projectives; return isMaximalNonCrossingCollection(k, n, union(list, collectionOfProjectives(k,n))); end
    return isMaximalNonCrossingCollection(k, n, list);
end

function isMaximalNonCrossingCollection(col::LabelCollection)
    return isMaximalNonCrossingCollection(col.k, col.n, col.collection)
end

function rotateSet(n, rotate, set)
    return sort(map(x->((x + rotate) % n == 0 ? n : (x+rotate) % n),set))
end

function rotateCollection(n, rotate, collection)
    return map(y->rotateSet(n,rotate,y),collection)
end

function rotate!(col::LabelCollection, rotate)
    col.collection = rotateCollection(col.n, rotate, col.collection)
end

"""
    Checks if a collection is symmetric.
"""
function isSymmetricCollection(k,n, collection)
    isSymmetricCollection(k, n, collection, false)
end

function isSymmetricCollection(col::LabelCollection)
    isSymmetricCollection(col.k, col.n, col.collection, false)
end

"""
    Checks if a collection is symmetric.
    Notice that this necessarily does not mean that the associated Postnikov diagram is symmetric.
"""
function isSymmetricCollection(k,n, collection, isSorted)
    rotatedCollection = rotateCollection(n, k, collection);
    if isSorted; return collection == sort(rotatedCollection); end
    return sort(collection) == sort(rotatedCollection)
end

"""
    Checks if two collections of labels col1 and col2 are equivalent up to rotation.
"""
function isEquivalentToCollectionUpToRotation(n, col1, col2)
    col1 = sort(col1)
    for i in 0:n
        c2 = sort(rotateCollection(n, i, col2))
        if col1 == c2; return true; end
    end
    return false
end

function isEquivalentToCollectionUpToRotation(col1::LabelCollection, col2::LabelCollection)
    if col1.n != col2.n && col1.k != col2.k
        return false
    end
    return isEquivalentToCollectionUpToRotation(col1.n, col1.collection, col2.collection)
end

"""
"""
function isEquivalentToCollectionUpToMirrorAndRotation(n, col1, col2)
    col1 = sort(col1)
    col1Mirr = sort(map(reverse, col1))
    for i in 0:n
        c2 = sort(rotateCollection(n, i, col2))
        if col1 == c2; return true; end
        if col1Mirr == c2; return true; end
    end
    return false
end

function isEquivalentToCollectionUpToMirrorAndRotation(col1::LabelCollection, col2::LabelCollection)
    if col1.n != col2.n && col1.k != col2.k
        return false
    end
    return isEquivalentToCollectionUpToRotation(col1.n, col1.collection, col2.collection)
end


"""
    
"""
function collectionOfProjectives(k,n)
    list = []
    for i in 1:n
        push!(list, sort(map(x-> x%n == 0 ? n : x % n , i:i+k-1)))
    end
    return list
end


#///////////////
# CReate Postnikov diagram

"""
    Generates a white clique with labels contining L.
    Make sure |L| = k - 1. With k being length of labels in collection
"""
function whiteClique(L, collection)
    list = []
    for label in collection
        if issubset(L, label)
            push!(list, label);
        end
    end
    return sort(list, by= x -> setdiff(x,L)[1]);
end


"""
    Generates a black clique with labels contained in L.
    Make sure |L| = k + 1. With k being length of labels in collection
"""
function blackClique(L, collection)
    list = []
    for label in collection
        if issubset(label, L)
            push!(list, label);
        end
    end
    return sort(list, by= x -> setdiff(L,x)[1]);
end

"""
    Returns a list of white cliques with length > 1, for a given noncrossing collections.
"""
function whiteCliques(collection)
    cliques = []
    labelsDone = []
    k = length(collection[1])
    for label in collection
        for l in combinations(label, k-1)
            if l in labelsDone continue end
            push!(labelsDone, l)
            wClique = whiteClique(l, collection)
            if length(wClique) > 2
                push!(cliques, wClique);
            end
        end
    end
    return cliques
end

function whiteCliques(col::LabelCollection)
    return whiteCliques(col.collection)
end

"""
    Returns a list of black cliques with length > 1, for a given noncrossing collections.
"""
function blackCliques(n, collection)
    cliques = []
    labelsDone = []
    for label in collection
        for i in 1:n
            if i in label; continue; end
            cliqueLabel = sort(union(label, [i]))
            if cliqueLabel in labelsDone
                continue
            end
            push!(labelsDone, cliqueLabel)
            bClique = blackClique(cliqueLabel, collection)
            if length(bClique) > 2
                push!(cliques, bClique)
            end
        end
    end
    return cliques;
end

function blackCliques(col::LabelCollection)
    return blackCliques(col.n, col.collection)
end

# What is this
# Path around clique 
# ex. [1,2,3,4] -> [[1,2], [2,3], [3,4], [4,1]]
# ex. [5,3,2,6] -> [[5,3], [3,2], [2,6], [6,5]]
function cliqueBoundary(clique, sortBy)
    bounds = []
    for i in 1:length(clique)
        j = (i == length(clique) ? 1 : i + 1)
        if sortBy == false
            push!(bounds, [clique[i], clique[j]])
            continue
        end
        push!(bounds, sort([clique[i], clique[j]], by=sortBy))
    end
    return bounds
end

function cliqueBoundary(clique)
    cliqueBoundary(clique, x->x);
end



"""
    findMaxinalNonCrossingCollections for not saving in file.
"""
function findMaxinalNonCrossingCollections(k, n)
    findMaxinalNonCrossingCollections(k, n, false);
end

"""
    Finding all maximal (k,n)-non-crossing collection of labels.
    NOTICE: This is done in a brute force way, making it very inefficient.
"""
function findMaxinalNonCrossingCollections(k, n, file; symmetric = true)
    if file != false
        open(file, "a") do io
            write(io, "\n")
            write(io, "*** Maximal (" * string(k) * "," * string(n) * ")-Non-Crossing collection of labels ***\n");
        end;
    end;
    list = []
    nonprojs = setdiff(collect(combinations(1:n,k)), collectionOfProjectives(k, n))
    combs = combinations(nonprojs,k*(n-k)-n+1)
    # Runs through all the combinations and check whether maximal
    println(length(combs))
    for m in combs
        if isMaximalNonCrossingCollection(k, n, m, true)
            if symmetric && !isSymmetricCollection(k, n, m, false)
                continue
            end

            isin=false
            for l in list
                if isEquivalentToCollectionUpToRotation(n, m, l); isin=true; break; end
            end
            if !isin;
                if file != false
                    open(file, "a") do io
                        write(io, string(m));
                        write(io, "\n");
                    end;
                end;
                push!(list, m); 
            end
        end
    end
    return list
end

function findMaxinalNonCrossingCollections2(k, n, file; symmetric = true)
    function next!(t,n,lastIndex)
        if t[lastIndex] == n 
            if lastIndex == 1
                t[lastIndex] = -1
                return
            end
            next!(t,n-1,lastIndex-1)
            t[lastIndex] = t[lastIndex-1]+1
        else
            t[lastIndex] = t[lastIndex] + 1
        end
    end

    if file != false
        open(file, "a") do io
            write(io, "\n")
            write(io, "*** Maximal (" * string(k) * "," * string(n) * ")-Non-Crossing collection of labels ***\n");
        end;
    end;
    list = []
    nonprojs = setdiff(collect(combinations(1:n,k)), collectionOfProjectives(k, n))
    combs = combinations(nonprojs,k*(n-k)-n+1)
    
    track = [i for i in 1:k]
    potentialCollection = []

    l = Int(n/k)
    

    function recF(tr, potential, add)
        for a in potential
            for b in add
                if !isNonCrossing(a,b)
                    return 
                end
            end
        end
        potential = union(potential, add)
        if length(potential) == k*(n-k)-n+1 # To make sure it is (k,n) - postnikov (otherwise can be other)
            if isMaximalNonCrossingCollection(k,n,union(potential,collectionOfProjectives(k, n)))
                if isSymmetricCollection(k, n, union(potential,collectionOfProjectives(k, n)), false)
                    println("SYM")
                    return;
                else
                    println("THERE IS A NON SYM")
                end
                #println("HERE")
                #println(potential)
                if file != false
                    isSymmetricCollection(k, n, union(potential,collectionOfProjectives(k, n)), false)
                    open(file, "a") do io
                        write(io, string(union(potential,collectionOfProjectives(k, n))));
                        write(io, "\n")
                    end;
                end;
                push!(list, potential);
                #return
            end
            return
        end
        t = copy(tr)
        while true
            next!(t,n,k)
            if t[begin] == -1
                return
            end
            if !(t in potential)
                if symmetric
                    recF(t,potential, [rotateSet(n,j*k, t) for j in 0:l-1])
                else
                    recF(t,potential, [t])
                end
            end
        end
    end

    recF(track, [], [])
end

