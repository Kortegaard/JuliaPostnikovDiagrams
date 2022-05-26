include("Postnikov.jl");
include("postnikovQuiver.jl");
include("maxcollections/maxcollection.jl")


k = 3
n = 9
cl = union(upToEquiv(k,n, col39))
println("There are " * string(length(cl)) * " diagram")

for i in 1:1
    println("Drawing number: " * string(i) * " / " * string(length(cl)))
    drawPostnikovDiagram(k,n,cl[i],
        filename                    = "test/quiv_"*string(k)*"_"*string(n)*"_" * string(i) * ".tex",
        showPostnikovQuiver         = false,
        showPlabicGraph             = true, 
        showPostnikovDiagram        = true, 
        showPostnikovDiamgramArrows = true,
    );
end


