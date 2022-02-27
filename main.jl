

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
            if !isNonCrossing(list[j],list[i]); return false; end
        end
    end
    return true
end

"""
    Checking if a collection of Noncorssing labels are maximal
"""
function isMaximalNonCrossingCollection(k, n, list)
    if length(list) == k*(n-k)+1; return isNonCrossingCollection(list); end;
    return false;
end

"""
    Checking if a collection of Noncorssing labels are maximal
"""
function isMaximalNonCrossingCollection(k, n, list, exclude_projectives)
    if exclude_projectives && length(list) == k*(n-k)-n+1; return isNonCrossingCollection(list); println("Checking..."); end;
    return isMaximalNonCrossingCollection(k, n, list);
end

function rotateCollection(n, rotate, collection)
    return map(y->sort(map(x->((x + rotate) % n == 0 ? n : (x+rotate) % n),y)),collection)
end

"""
    Checks if a collection is symmetric.
"""
function isSymmetricCollection(k,n, collection)
    isSymmetricCollection(k, n, collection, false)
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

"""
    
"""
function collectionOfProjections(k,n)
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

using Combinatorics

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
function findMaxinalNonCrossingCollections(k, n, file)
    if file != false
        open(file, "a") do io
            write(io, "\n")
            write(io, "*** Maximal (" * string(k) * "," * string(n) * ")-Non-Crossing collection of labels ***\n");
        end;
    end;
    list = []
    nonprojs = setdiff(collect(combinations(1:n,k)), collectionOfProjections(k, n))
    combs = combinations(nonprojs,k*(n-k)-n+1)
    for m in combs
        if isMaximalNonCrossingCollection(k, n, m, true) && isSymmetricCollection(k, n, m, false)
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


#println(findMaxinalNonCrossingCollections(2,6, "out.txt"));

#col48 = union([[1, 2, 3, 5], [1, 2, 3, 6], [1, 2, 3, 7], [1, 2, 4, 5], [1, 2, 5, 6], [1, 5, 6, 7], [1, 5, 6, 8], [2, 5, 6, 7], [3, 5, 6, 7]], collectionOfProjections(4,8))
include("maxcollection.jl")
#col48_2 = union([[ 1, 2, 3, 5], [1, 2, 3, 6], [1, 2, 3, 7], [1, 2, 5, 6], [1, 3, 4, 5], [1, 5, 6, 7], [1, 5, 7, 8], [2, 5, 6, 7], [3, 5, 6, 7]],collectionOfProjections(4,8))
col48_2 = union(col48[40],collectionOfProjections(4,8))
#col48_2 = union(col48[19],collectionOfProjections(4,8))


l1 = whiteCliques(col48_2)
l2 = blackCliques(8, col48_2)
#open("Maximal ()")

println(col48_2)
println()
println("** WHITE CLIQUES");
for m in l1
    println(m)
end
println();
println("** BLACK CLIQUES");
for m in l2
    println(m)
end

notins = setdiff(collect(combinations(1:n,k)), col48_2)

#println(isMaximalNonCrossingCollection(4,8, col48_2))
#col = union([[1, 2, 4], [1, 2, 5], [1, 4, 5], [2, 4, 5]],collectionOfProjections(3,6))
#println(whiteCliques(col))j
#println(blackCliques(6, col))

#println(whiteClique([1,2], union([[1, 2, 4], [1, 2, 5], [1, 4, 5], [2, 4, 5]],collectionOfProjections(3,6))))