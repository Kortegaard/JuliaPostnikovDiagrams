include("/Users/ank/master/code/julia_quiver/src/main.jl")
include("/Users/ank/master/code/julia_quiver/src/plot.jl")
using LinearAlgebra

function quiverFromCollection(k,n,collection)::Quiver
    q = Quiver(map(x->Vertex(join(x)),collection))#, ["a",1,2], ["b",3,2], ["c",2,4])

    wc = whiteCliques(collection)
    bc = blackCliques(n,collection)
    wc_labels = map(y->map(x->join(x),y),wc)
    bc_labels = map(y->map(x->join(x),y),bc)
    for clique in wc_labels
        for i in 1:length(clique)
            j = (i == length(clique) ? 1 : i + 1)
            if !exists_arrow(q, Vertex(clique[i]), Vertex(clique[j]))
                add_arrow!(q, Vertex(clique[i]), Vertex(clique[j]))
            end
        end
    end
    for clique in bc_labels
        for i in 1:length(clique)
            j = (i == length(clique) ? 1 : i + 1)
            if !exists_arrow(q, Vertex(clique[j]), Vertex(clique[i]))
                add_arrow!(q, Vertex(clique[j]), Vertex(clique[i]))
            end
        end
    end

    #using Plots
    set_random_positions!(q);
    t = 0
    for mm in collectionOfProjectives(k,n)
        v = find_vertex(q,Vertex(join(sort(mm))))
        v.data["springFrozen"] = true;
        v.data["position"] = [4*cos(2*pi*t/n),4*sin(2*pi*t/n)]
        t = t+1
    end

    for i in 1:1000
        spring_step(q,0.1,1.0,1.0)
    end

    normalize_quiver!(q)

    return q
end

function mean_vector(vects::Vector{Vector{Float64}})
    return sum(vects)/length(vects)
end

function constructCliqueQuiver(k,n, collection, collectionQuiver)
    wc = whiteCliques(collection)
    bc = blackCliques(n,collection)
    cliqueQuiver = Quiver(union([Vertex("w"*string(i)) for i in 1:length(wc)], [Vertex("b"*string(i)) for i in 1:length(bc)]))

    for i in 1:length(wc)
        v = find_vertex(cliqueQuiver, Vertex("w"*string(i)))
        v.data = Dict{String, Any}()
        v.data["clique"] = [find_vertex(collectionQuiver, Vertex(join(x))) for x in wc[i]]
        v.data["cliqueLabels"] = [join(x) for x in wc[i]]
        v.data["position"] = mean_vector([find_vertex(collectionQuiver,Vertex(join(x))).data["position"] for x in wc[i]])
        v.data["color"] = "white"
    end
    for i in 1:length(bc)
        v = find_vertex(cliqueQuiver, Vertex("b"*string(i)))
        v.data = Dict{String, Any}()
        v.data["clique"] = [find_vertex(collectionQuiver, Vertex(join(x))) for x in bc[i]]
        v.data["cliqueLabels"] = [join(x) for x in bc[i]]
        v.data["position"] = mean_vector([find_vertex(collectionQuiver,Vertex(join(x))).data["position"] for x in bc[i]])
        v.data["color"] = "black"
    end

    for c1 in vertices(cliqueQuiver)
        for c2 in vertices(cliqueQuiver)
            if c1 == c2; continue; end
            int = intersect(cliqueBoundary(c1.data["clique"], x->x.label), cliqueBoundary(c2.data["clique"], x->x.label))
            if length(int) == 1
                add_arrow!(cliqueQuiver,c1,c2)
            end
        end
    end

    for i in 1:n
        v = Vertex(string(i))
        v.data = Dict{String, Any}()
        offset = 0
        if n%2 == 0
            offset = n+3 # n-1
        else
            offset = n+2
        end
        v.data["position"] = [cos(2*pi*(2*i + (offset))/(2*n)), sin(2*pi*(2*i + (offset))/(2*n))]
        add_vertex!(cliqueQuiver, v)
    end

    for b in cliqueBoundary(collectionOfProjectives(k,n), false)
        blab = sort(map(join, sort(b)))
        v = find_vertex(cliqueQuiver, Vertex(string(setdiff(b[2], b[1])[1])))
        for vv in vertices(cliqueQuiver)
            if !haskey(vv.data, "cliqueLabels") continue end
            if blab in cliqueBoundary(vv.data["cliqueLabels"])
                add_arrow!(cliqueQuiver, vv, v)
            end
        end
    end
    return cliqueQuiver
end 

import Base:angle
function angle(a, b)
    aa = [-a[2], a[1]]
    sgn = sign(dot(aa, b))
    sgn = (sgn == 0 ? 1 : sgn)
    return sgn*acosd(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
end

#
function dirArrowFromToVertex(q::Quiver, f::Vertex, t::Vertex, dir::String)
    neighbours_verts = neighbours(q,f)
    setdiff!(neighbours_verts,[t])
    if length(neighbours_verts) == 0 
        return nothing
    end
    diff = t.data["position"] - f.data["position"]
    neighbours_angles = [angle(diff, neighbours_verts[i].data["position"] - f.data["position"] ) for i in 1:length(neighbours_verts)]
    neighbours_angles =  map(x -> (x < 0 ? x + 360 : x), neighbours_angles)
    #println(neighbours_angles)

    if dir == "right"
        return neighbours_verts[argmax(neighbours_angles)]
    end
    if dir == "left"
        return neighbours_verts[argmin(neighbours_angles)]
    end
    return nothing
end




##################3

using Interpolations

function postnikovDiagramDrawingData(cliqueQuiver, n)
    lns = []
    for i in 1:n
        v = find_vertex(cliqueQuiver, Vertex(string(i)))
        neighs = neighbours(cliqueQuiver, v)
        if length(neighs) == 0; continue; end
        n = neighs[1]

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
            

            nn = dirArrowFromToVertex(cliqueQuiver, n, v, nn_dir)
            if nn == nothing 
                break
            end
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







########## PLOTTING ###########

function plot_quiver(qq::Quiver; directed=true, vertex_color="black", linecolor="black", arrow_offset = 0.1)
    #Arrows
    for arr in arrows(qq)
            offset = arrow_offset
            offset2 = (1-offset)^2
            dir = arr.termination.data["position"]-arr.start.data["position"]
            if directed
                arrows!(map(x->[x],arr.start.data["position"]+dir*offset)..., map(x->[x], arr.termination.data["position"]-arr.start.data["position"])..., lengthscale=offset2, color=linecolor, linewidth=1.5)
            else
                lines!( map(x->[i for i in x],collect(zip(arr.start.data["position"], arr.termination.data["position"]) ))..., color=linecolor)
            end
    end

    #Points
    scatter!(map(x -> [i for i in x],collect(zip(map(x->x.data["position"],vertices(qq))...)))..., legend = false, color=vertex_color, markersize=10)

end

function drawPostnikovDiagram(k,n,maximalNonCrossingCollection; fig = nothing, showPlabicGraph = false, showPostnikovDiagram = true, showPostnikovDiamgramArrows = true, showPostnikovQuiver=false, saveAs=false, drawOuterCirle=true)

    postnikovQuiver = quiverFromCollection(k,n, maximalNonCrossingCollection);
    cliqueQuiver = constructCliqueQuiver(k,n, maximalNonCrossingCollection, postnikovQuiver)
    m, arrs = postnikovDiagramDrawingData(cliqueQuiver, n)

    if fig == nothing
        fig = Figure(resolution = (800 , 800), backgroundcolor = :transparent,figure_padding = 1);
        ax = Axis(fig[1, 1], backgroundcolor = :transparent)
        hidedecorations!(ax)
        hidespines!(ax)
        ax.xticksvisible = false
        ax.yticksvisible = false
        ax.backgroundcolor = :transparent
    #else
        #ax = Axis(fig[1, 1], backgroundcolor = :transparent)
    end

    # Position Postnikov vertices in the middle of cliques
    #for v in vertices(postnikovQuiver)
    #    if !haskey(v.data, "springFrozen") || v.data["springFrozen"] == false
    #        v_clique_vertices = filter(x-> haskey(x.data, "clique") && v in x.data["clique"], vertices(cliqueQuiver))
    #        mean_vec = mean_vector(map(x->x.data["position"], v_clique_vertices))
    #        v.data["position"] = mean_vec
    #    end
    #end

#map(x->x.data["position"], m)

    if drawOuterCirle
        sides = 100
        rad = 1
        circle_xs = [rad*cos(2*pi*i/sides) for i in 0:sides]
        circle_ys = [rad*sin(2*pi*i/sides) for i in 0:sides]
        lines!(circle_xs, circle_ys, linewidth=2, linestyle = :dash)
    end


    #println("HE")

    if showPlabicGraph
        plot_quiver(cliqueQuiver, directed=false, vertex_color=:red);
    end
    if showPostnikovQuiver
        plot_quiver(postnikovQuiver, directed=true);
    end

    if showPostnikovDiagram
        for mm in m
            lines!( mm[1], mm[2], color="green")
        end
    end
    if showPostnikovDiagram && showPostnikovDiamgramArrows
        for arr in arrs
            scatter!([arr[1][1]], [arr[1][2]], color="green", marker='▲', rotations=[arr[2]])
        end
    end
    if saveAs
        save(saveAs, f)
    end

    return fig

end

