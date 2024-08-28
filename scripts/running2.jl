include("../src/collections.jl");
include("../src/postnikovQuiver.jl");
include("../src/mutation.jl");
include("../maxcollections/maxcollection.jl")
include("../maxcollections/out312.jl")



k = 3
n = 9
c39 = union(upToEquiv(k,n, col39))
NC = LabelCollection(3,9,c39[1])
NCmy = LabelCollection(3,9,[[1,2,3],[2,3,4],[3,4,5],[4,5,6],[5,6,7],[6,7,8],[7,8,9],[1,8,9],[1,2,9],
                            [2,8,9],[2,3,9],[3,8,9],[6,8,9],[5,6,8],[5,6,9],[3,6,9],[2,3,6],[3,5,6],[2,3,5]])

mut1 = mutationAtLabel(NCmy,[2,8,9])[1]
mut2 = mutationAtLabel(mut1,[5,6,8])[1]
mut3 = mutationAtLabel(mut2,[2,3,5])[1]
cl = mutations(NC)
# println(whiteCliques(col36[1]))
# println(blackCliques(6, col36[1]))

#cl = union(upToEquiv(k,n, col36))
#println("there are " * string(length(cl)) * " diagram")

drawPostnikovDiagram(k,n,NCmy.collection,
    filename                    = "test/start.tex",
    showPostnikovQuiver         = false,
    showPlabicGraph             = false, 
    showPostnikovDiagram        = true, 
    showPostnikovDiamgramArrows = true,
);
drawPostnikovDiagram(k,n,mut1.collection,
    filename                    = "test/mut1-239.tex",
    showPostnikovQuiver         = false,
    showPlabicGraph             = false, 
    showPostnikovDiagram        = true, 
    showPostnikovDiamgramArrows = true,
);
drawPostnikovDiagram(k,n,mut2.collection,
    filename                    = "test/mut2-568.tex",
    showPostnikovQuiver         = false,
    showPlabicGraph             = false, 
    showPostnikovDiagram        = true, 
    showPostnikovDiamgramArrows = true,
);
drawPostnikovDiagram(k,n,mut3.collection,
    filename                    = "test/mut3-235.tex",
    showPostnikovQuiver         = false,
    showPlabicGraph             = false, 
    showPostnikovDiagram        = true, 
    showPostnikovDiamgramArrows = true,
);




#for i in 1:length(cl)
#    println("drawing number: " * string(i) * " / " * string(length(cl)))
#    drawPostnikovDiagram(k,n,cl[i].collection,
#        filename                    = "test/mpd_"*string(k)*"_"*string(n)*"_" * string(i) * ".tex",
#        showPostnikovQuiver         = false,
#        showPlabicGraph             = false, 
#        showPostnikovDiagram        = true, 
#        showPostnikovDiamgramArrows = true,
#    );
#end


