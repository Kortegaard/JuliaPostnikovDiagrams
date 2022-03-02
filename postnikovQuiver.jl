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
        v.data["position"] = [cos(2*pi*(2*i + (n-1))/(2*n)), sin(2*pi*(2*i + (n-1))/(2*n))]
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


function plot_quiver(fig::Figure, qq::Quiver)
    plot_quiver(fig, qq, true)
end

function plot_quiver(fig::Figure, qq::Quiver, directed::Bool)
    ax = Axis(fig[1, 1], backgroundcolor = :transparent)
    hidedecorations!(ax)
    hidespines!(ax)

    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.backgroundcolor = :transparent


    #cirlce
    sides = 100
    rad = 1
    circle_xs = [rad*cos(2*pi*i/sides) for i in 0:sides]
    circle_ys = [rad*sin(2*pi*i/sides) for i in 0:sides]
    lines!(circle_xs, circle_ys, linewidth=2, linestyle = :dash)
    ###

    #Arrows
    for arr in arrows(qq)
            offset = 0.1
            offset2 = (1-offset)^2
            dir = arr.termination.data["position"]-arr.start.data["position"]
            if directed
                arrows!(map(x->[x],arr.start.data["position"]+dir*offset)..., map(x->[x], arr.termination.data["position"]-arr.start.data["position"])..., lengthscale=offset2)
            else
                lines!( map(x->[i for i in x],collect(zip(arr.start.data["position"], arr.termination.data["position"]) ))..., color="black")
            end
            #println( map(x->[i for i in x],collect(zip(arr.start.data["position"], arr.termination.data["position"]) )) )

    end

    #Points
    scatter!(map(x -> [i for i in x],collect(zip(map(x->x.data["position"],vertices(qq))...)))..., legend = false)

end

function save_quiver_plot(q::Quiver, fileName::String)
    f = Figure(resolution = (800 , 800), backgroundcolor = :transparent,figure_padding = 1)
    plot_quiver(f, q)
    save(fileName, f)
end















#open("test.gn", "w") do io
#    write(io, "set xrange [-5:5]\n")
#    write(io, "set yrange [-5:5]\n")
#    write(io, "set style arrow 1 head filled size screen 0.02,15,15 ls 1\n")
#    write(io, "plot NaN t ''\n")
#    for arr in arrows(q)
#        write(io, "set arrow from ") 
#        write(io, string(arr.start.data["position"][1]),",",string(arr.start.data["position"][2]))
#        write(io, " to ") 
#        write(io, string(arr.termination.data["position"][1]),",",string(arr.termination.data["position"][2]))
#        write(io, " as 1\n")
#    end
#    write(io, "plot '-' w points lc rgb \"black\" pt 7\n")
#    for v in vertices(q)
#        write(io, string(v.data["position"][1]), " ", string(v.data["position"][2], "\n"))
#    end
#    write(io, "EOF\n")
    #write(io, "plot NaN t ''\n")
#
#end

