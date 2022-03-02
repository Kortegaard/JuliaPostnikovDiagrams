include("postnikov.jl");
include("maxcollection.jl")



#liss = []
#for m in col48
#    if isMaximalNonCrossingCollection(4,8,union(m,collectionOfProjectives(4,8)))
#        push!(liss, union(m,collectionOfProjectives(4,8)))
#    end
#end
#
#liss2 = []
#for m in col39
#    if isMaximalNonCrossingCollection(3,9,union(m,collectionOfProjectives(3,9)))
#        push!(liss2, union(m,collectionOfProjectives(3,9)))
#    end
#end

####################3

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


n = 10
k = 5
function findMaxinalNonCrossingCollections(k, n, file)
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
        if length(potential) == k*(n-k)-n+1
            if isMaximalNonCrossingCollection(k,n,potential)
                #println(potential)
                push!(list, potential);
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
                recF(t,potential, [rotateSet(n,j*k, t) for j in 0:l-1])
            end
        end
    end

    recF(track, [], [])

end

q = quiverFromCollection(5,10,col510_p[1])
cQ = constructCliqueQuiver(5,10, col510_p[1], q)


using Interpolations
function PostnikovDiagramLines(cliqueQuiver, n)
    lns = []
    for i in 1:n
        v = find_vertex(cliqueQuiver, Vertex(string(i)))
        n = neighbours(cliqueQuiver, v)[1]

        points =  [v.data["position"]]
        # Replace 100 by a more sensable choise ()
        for _ in 1:100
            nn_dir = ""
            if haskey(n.data, "color") && n.data["color"] == "black"
                nn_dir = "right"
            elseif haskey(n.data, "color") && n.data["color"] == "white"
                nn_dir = "left"
            else
                break
            end
            

            nn = dirArrowFromToVertex(cQ, n, v, nn_dir)
            v = n
            n = nn
            if haskey(n.data, "color")
                push!(points, mean_vector([n.data["position"], v.data["position"]]))
            else
                push!(points, n.data["position"])
            end
        end

        # Making the lines smoother using cubic interpolation (Interpolation library)
        A = hcat(map(x->[i for i in x], collect(zip(points...)))...)
        t = 0:(1/(length(points)-1)):1
        itp = Interpolations.scale(interpolate(A, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)
        tfine = 0:.01:1
        xs, ys = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]

        push!(lns, [xs,ys,length(points)])
    end

    arrs = []
    for mm in lns
        
        for j in 0:mm[3]-2
        #for j in 0:0
            index = ((2*j+1)*(length(mm[1]))÷(2*(mm[3]-1)))+1
            point = [mm[1][index], mm[2][index]]
            next_point = [mm[1][index+2], mm[2][index+2]]
            dir = (next_point - point)
            agl = angle([0.0,1.0], dir)*pi /180
            
            push!(arrs, [point, agl]);
        end
    end

    return lns, arrs
end

f = Figure(resolution = (800 , 800), backgroundcolor = :transparent,figure_padding = 1);
ax = Axis(f[1, 1], backgroundcolor = :transparent)

hidedecorations!(ax)
hidespines!(ax)
ax.xticksvisible = false
ax.yticksvisible = false
ax.backgroundcolor = :transparent

plot_quiver(f, q);
plot_quiver(f, cQ,false);
m, arrs = PostnikovDiagramLines(cQ, 10)
for mm in m
    lines!( mm[1], mm[2], color="green")
end
for arr in arrs
    scatter!([arr[1][1]], [arr[1][2]], color="green", marker='▲', rotations=[arr[2]])
end
save("fn4.pdf", f)

