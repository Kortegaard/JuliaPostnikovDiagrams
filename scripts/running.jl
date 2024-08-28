include("../src/collections.jl");
include("../src/postnikovQuiver.jl");
include("../src/mutation.jl");
include("../maxcollections/maxcollection.jl")
include("../maxcollections/out312.jl")


k = 3
n = 6

NC = LabelCollection(3,6,col36[1])
println(NC)
cl = mutations(NC)
# println(whiteCliques(col36[1]))
# println(blackCliques(6, col36[1]))

#cl = union(upToEquiv(k,n, col36))
#println("there are " * string(length(cl)) * " diagram")


for i in 1:length(cl)
    println("drawing number: " * string(i) * " / " * string(length(cl)))
    drawPostnikovDiagram(k,n,cl[i].collection,
        filename                    = "test/mpd_"*string(k)*"_"*string(n)*"_" * string(i) * ".tex",
        showPostnikovQuiver         = false,
        showPlabicGraph             = false, 
        showPostnikovDiagram        = true, 
        showPostnikovDiamgramArrows = true,
    );
end


