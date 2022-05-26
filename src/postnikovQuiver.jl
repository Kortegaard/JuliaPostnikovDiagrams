include("../julia_quiver/src/main.jl")
include("../julia_quiver/src/plot.jl")
using LinearAlgebra

mutable struct PostnikovDiagram
    # For a (k,n) - Postnikov diagram
    k::Int
    n::Int

    quiver::Quiver # The associated quiver to the postnikov diagram
    plabicGraph::Quiver

    whiteCliques::Vector{Vector{Arrow}}
    blackCliques::Vector{Vector{Arrow}}

    collection::Vector{Vector{Any}}
end


"""
    Constructing the Postnikov data given a maximal collection.
"""
function PostnikovDiagram(k,n,maxNonCrossColl)
    q = quiverFromCollection(k,n,maxNonCrossColl);
    cQ = constructCliqueQuiver(k, n, maxNonCrossColl, q);


    wc = whiteCliques(maxNonCrossColl)
    bc = blackCliques(n, maxNonCrossColl)

    wl = []
    for cc in wc
        l1cliq = []
        for i in 1:length(cc)
            j = (i+1 > length(cc) ? 1 : i+1)
            v1 = find_vertex(q, Vertex(join(cc[i])))
            v2 = find_vertex(q, Vertex(join(cc[j])))
            a = find_arrow(q, v1,v2)
            push!(l1cliq, a);
        end
        push!(wl, l1cliq);
    end

    bl = []
    for cc in bc
        l2cliq = []
        for i in 1:length(cc)
            j = (i-1 < 1 ? length(cc) : i-1)
            v1 = find_vertex(q, Vertex(join(cc[i])))
            v2 = find_vertex(q, Vertex(join(cc[j])))
            a = find_arrow(q, v1,v2)
            push!(l2cliq, a);
        end
        push!(bl, l2cliq);
    end

   return PostnikovDiagram(k,n, q, cQ, wl,bl, maxNonCrossColl)
end

function quiverFromCollection(k,n,collection)::Quiver
    q = Quiver(map(x->Vertex(join(x)),collection)) #, ["a",1,2], ["b",3,2], ["c",2,4])

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
        v.data["frozen"] = true;
        v.data["position"] = [4*cos(2*pi*t/n),4*sin(2*pi*t/n)]
        t = t+1
    end

    for i in 1:1000
        spring_step(q,0.1,0.2,1.0)
        #spring_step(q,0.1,1.0,1.0)
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
            offset = n+5 #n+3 # n-1
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
    scatter!(map(x -> [i for i in x],collect(zip(map(x->x.data["position"],vertices(qq))...)))..., color=vertex_color, markersize=10)

end


function tikz_plot_quiver(qq::Quiver; directed=true, vertex_color="black", linecolor="black", arrow_offset = 0.03, linewidth = 0.8, markerSize=1, draw_vertices=true ,kwargs...)
    #Arrows
    output = ""
    for arr in arrows(qq)
        if directed
            output *= tikzDrawArrow(arr.start.data["position"], arr.termination.data["position"], offset=arrow_offset, linewidth=linewidth)
        else
            output *= tikzDrawLine( map(x->[i for i in x],collect(zip(arr.start.data["position"], arr.termination.data["position"]) ))..., color=linecolor, linewidth=linewidth)
        end
    end

    #Points
    if draw_vertices
        output *= tikzPoints(map(x -> [i for i in x],collect(zip(map(x->x.data["position"],vertices(qq))...)))..., color=vertex_color, markSize=1.6)
    end
    return output

end

function tikzDrawCommand(;linewidth=0.05, color="black", dashed=false, dashPattern=(1.0,1.0), extra="")
    output = "\\draw["
    output = output * "line width="*string(linewidth)*"pt, "
    #if dashed
    #    output = output * "dashed, "
    #end
    if dashed
        output = output * "dash pattern=on "*string(dashPattern[1])*"pt off "*string(dashPattern[2])*"pt, "
    end
    output = output * string(color) * ","
    output *= extra
    output *= "] "
    return output
end

function tikzDrawCirle(center, radius; kwargs...)
    output = tikzDrawCommand(;kwargs...) * " ("*string(center[1])*","*string(center[2])*") circle ("*string(radius)*");"
    return output
end

function tikzDrawLine(xs, ys; kwargs...)
    output = tikzDrawCommand(;kwargs...) * " "
    points = [[xs[i],ys[i]] for i in 1:length(xs)]
    for i in 1:length(points)
        output = output * "("*string(points[i][1])*","*string(points[i][2])*")"
        if i != length(points)
            output = output * " -- "
        else
            output = output * ";\n"
        end
    end
    return output
end

function tikzDrawArrow(from, to; offset=0, kwargs...)
    output = tikzDrawCommand(extra="-{Stealth[length=6pt]},"; kwargs...) * " "
    vec = (to - from);
    vec = offset*(vec/norm(vec))
    from = from + vec;
    to = to - vec
    output *= "("*string(from[1])*","*string(from[2])*")" * " -- " * "("*string(to[1])*","*string(to[2])*");\n";
    return output
end

function tikzPoint(x, y; rotate=0, color="black", marker="*", markSize=0.1, kwargs...)
    output = ""
    output = output * "\\node[mark size="*string(markSize)*"pt,color="*color*",rotate="*string(rotate)*"] at ("*string(x)*","*string(y)*") {\\pgfuseplotmark{"*marker*"}}"
    output = output * ";\n"
    return output
end

function tikzPoints(xs, ys; kwargs...)
    output = ""
    for i in 1:length(xs)
        output *= tikzPoint(xs[i], ys[i]; kwargs...)
    end
    return output
end

function drawPostnikovDiagram(k,n,maximalNonCrossingCollection;filename="", fig = nothing, showPlabicGraph = false, showPostnikovDiagram = true, showPostnikovDiamgramArrows = true, showPostnikovQuiver=false, saveAs=false, drawOuterCirle=true)

    open(filename, "w") do file
        write(file, "\\documentclass[crop,tikz]{standalone}\n\\usetikzlibrary{plotmarks,arrows.meta}\n\\definecolor{mycolor}{rgb}{0.152941, 0.682353, 0.376471}\\begin{document}\n\\begin{tikzpicture}[x=150pt,y=150pt]\n") 
        postnikovQuiver = quiverFromCollection(k,n, maximalNonCrossingCollection);
        cliqueQuiver = constructCliqueQuiver(k,n, maximalNonCrossingCollection, postnikovQuiver)
        m, arrs = postnikovDiagramDrawingData(cliqueQuiver, n)

        if drawOuterCirle
            write(file, tikzDrawCirle((0,0), 1, color="blue", dashPattern=(3.0,3.0), dashed=true, linewidth=0.8))
        end

        if showPlabicGraph
            write(file,tikz_plot_quiver(cliqueQuiver, directed=false, vertex_color="red", linewidth=0.8, draw_vertices=false));
            out = ""
            for v in vertices(cliqueQuiver)
                if haskey(v.data, "color")
                    if getindex(v.data, "color") == "black"
                        out *= tikzPoint(v.data["position"][1],v.data["position"][2], color="black", markSize=1.6)
                    end
                    if getindex(v.data, "color") == "white"
                        out *= tikzPoint(v.data["position"][1],v.data["position"][2], color="red", markSize=1.6)
                    end
                else
                    out *= tikzPoint(v.data["position"][1],v.data["position"][2], color="blue", markSize=1.6)
                end
            end
            write(file, out);
        end

        if showPostnikovQuiver
            write(file,tikz_plot_quiver(postnikovQuiver, directed=true));
        end

        if showPostnikovDiagram
            for mm in m
                write(file,tikzDrawLine(mm[1],mm[2], linewidth=0.6, color="mycolor"))
            end
        end

        if showPostnikovDiagram && showPostnikovDiamgramArrows
            for arr in arrs
                write(file,tikzPoint(arr[1][1], arr[1][2], rotate=(arr[2]/pi * 180), marker="triangle*", markSize=1.6, color="mycolor"))
            end
        end

        if saveAs
            save(saveAs, f)
        end

        write(file, "\\end{tikzpicture}\n\\end{document}")
    end

    return fig
end
