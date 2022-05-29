include("postnikovQuiver.jl")

using GAP
GAP.evalstr("LoadPackage(\"QPA\");")

"""
    Given a quiver q, this function will return a string of a GAP QPA quiver.

    TODO: Change `includeFrozen`, to a filter function.
"""
function quiverToQPAString(q::Quiver; includeFrozen = true)
    i = 1
    for v in vertices(q)
        if !haskey(v.data, "name")
            v.data["name"] = string(i)
        end
        i = i+1
    end
    i = 1
    for a in arrows(q)
        if !haskey(a.data, "name")
            a.data["name"] = "a" * string(i)
        end
        i = i+1
    end

    function vertexIsFrozen(v::Vertex)::Bool
        return haskey(v.data, "frozen") && v.data["frozen"]==true
    end
    vertsString = string([v.data["name"] for v in vertices(q) if (includeFrozen || !vertexIsFrozen(v))])
    arrString = string([[a.start.data["name"], a.termination.data["name"], a.data["name"]] for a in arrows(q) if (includeFrozen || (!vertexIsFrozen(a.start) && !vertexIsFrozen(a.termination)))])
    
    gapQuiverStr = "Quiver(" * vertsString * ", " * arrString * ");"
    return gapQuiverStr
end


"""
    Given a Postnikov diagram pd, this function will generate the ideal of which comes from the potential
    of the frozen quiver associated to Postnikov diagram pd
"""
function qpaPostnikovQuiverIdeal(pd::PostnikovDiagram; includeFrozen = true, pre = "Q.")
    idealList = []
    for arr in arrows(pd.quiver)
        if !includeFrozen
            if get(arr.start.data, "frozen", false) || get(arr.termination.data, "frozen", false)
                continue
            end
        end

        wcEl = get(filter(x -> arr in x, pd.whiteCliques), 1,nothing)
        if any(map(arr->get(arr.start.data, "frozen", false) || get(arr.termination.data, "frozen", false),wcEl))
            wcEl = nothing
        end
        bcEl = get(filter(x -> arr in x, pd.blackCliques), 1,nothing)
        if any(map(arr->get(arr.start.data, "frozen", false) || get(arr.termination.data, "frozen", false),bcEl))
            bcEl = nothing
        end

        wStr = ""
        if wcEl != nothing 
            wff = findfirst(z->z==arr, wcEl)
            wrestOfCirc = [wcEl[((wff+i) %length(wcEl)+ 1)]  for i in 0:length(wcEl)-2]
            wStr  = join((map(x->pre*x.data["name"], wrestOfCirc)), "*")
        end

        bStr = ""
        if bcEl != nothing 
            bff = findfirst(z->z==arr, bcEl)
            brestOfCirc = [bcEl[((bff+i) %length(bcEl)+ 1)]  for i in 0:length(bcEl)-2]
            bStr  = "-"*join(reverse(map(x->pre*x.data["name"], brestOfCirc)), "*")
        end

        str = wStr*bStr
        push!(idealList, str)
    end
    return "[" * join(idealList, ", ") * "]"
end

#for myCol in upToEquiv(3,9,col39)# col412
#    mpd = PostnikovDiagram(3,9,myCol);
#    Q = quiverToQPAString(mpd.quiver, includeFrozen=false)
#    sss = qpaPostnikovQuiverIdeal(mpd, includeFrozen=false, pre="kQ.")
#
#    GAP.evalstr("Q:="*Q*";")
#    GAP.evalstr("kQ := PathAlgebra(GF(3), Q);")
#    GAP.evalstr("I := Ideal(kQ,"*sss*");")
#    GAP.evalstr("IsAdmissibleIdeal(I);")
#    GAP.evalstr("A := kQ/I;")
#
#    println(GAP.evalstr("Determinant(CartanMatrix(A));"))
#    println(GAP.evalstr("IsSymmetricAlgebra(A);"))
#    println("")
#end