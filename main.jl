include("postnikov.jl");
include("postnikovQuiver.jl");
include("maxcollections/maxcollection.jl")
#include("maxcollections/out412.jl")
#include("maxcollections/out312.jl")
#include("maxcollections/out510.jl")

#function drawFourDiags(k,n,coll, fileName)
    ##fig = Figure(resolution = (1600, 1600), backgroundcolor = :transparent,figure_padding = 1);
    #fig = Figure(resolution = (800,800), backgroundcolor = :transparent,figure_padding = 1);

    #ax = Axis(fig[1, 1], backgroundcolor = :transparent)
    #hidedecorations!(ax)
    #hidespines!(ax)
    #ax.xticksvisible = false
    #ax.yticksvisible = false
    #ax.backgroundcolor = :transparent
    #poshDrag = drawPostnikovDiagram(k,n,coll,
        #showPostnikovQuiver         = false,
        #showPlabicGraph             = false, 
        #showPostnikovDiagram        = true, 
        #showPostnikovDiamgramArrows = true,
        #fig = fig
    #);

    ##ax2 = Axis(fig[1, 2], backgroundcolor = :transparent)
    ##hidedecorations!(ax2)
    ##hidespines!(ax2)
    ##ax2.xticksvisible = false
    ##ax2.yticksvisible = false
    ##ax2.backgroundcolor = :transparent
    ##poshDrag = drawPostnikovDiagram(k,n,coll,
        ##showPostnikovQuiver         = true,
        ##showPlabicGraph             = false, 
        ##showPostnikovDiagram        = false, 
        ##showPostnikovDiamgramArrows = true,
        ##fig = fig
    ##);
    

    ##ax3 = Axis(fig[2, 1], backgroundcolor = :transparent)
    ##hidedecorations!(ax3)
    ##hidespines!(ax3)
    ##ax3.xticksvisible = false
    ##ax3.yticksvisible = false
    ##ax3.backgroundcolor = :transparent
    ##poshDrag = drawPostnikovDiagram(k,n,coll,
        ##showPostnikovQuiver         = true,
        ##showPlabicGraph             = false, 
        ##showPostnikovDiagram        = true, 
        ##showPostnikovDiamgramArrows = true,
        ##fig = fig
    ##);


    ##ax4 = Axis(fig[2, 2], backgroundcolor = :transparent)
    ##hidedecorations!(ax4)
    ##hidespines!(ax4)
    ##ax4.xticksvisible = false
    ##ax4.yticksvisible = false
    ##ax4.backgroundcolor = :transparent
    ##poshDrag = drawPostnikovDiagram(k,n,coll,
        ##showPostnikovQuiver         = true,
        ##showPlabicGraph             = true, 
        ##showPostnikovDiagram        = false, 
        ##showPostnikovDiamgramArrows = true,
        ##fig = fig
    ##);


    #save(fileName, fig);
#end


k = 3
n = 9
##findMaxinalNonCrossingCollections2(3,9,"./out39_2.txt")
cl = union(upToEquiv(k,n, col39))
##cl = co
##
#println("BEGINNING")
#findMaxinalNonCrossingCollections2(3,9,"./out.txt", symmetric=false)
##findMaxinalNonCrossingCollections(3,6,"./36file.txt");
#println("There are " * string(length(cl)) * " diagram")
for i in 1:1
    println("Drawing number: " * string(i) * " / " * string(length(cl)))
    drawPostnikovDiagram(k,n,cl[i],
        filename                    = "test/quiv_"*string(k)*"_"*string(n)*"_" * string(i) * ".tex",
        showPostnikovQuiver         = false,
        showPlabicGraph             = true, 
        showPostnikovDiagram        = true, 
        showPostnikovDiamgramArrows = true,
    );
   # drawPostnikovDiagram(k,n,cl[i],
   #     filename                    = "test/quiv_"*string(k)*"_"*string(n)*"_" * string(i) * "_2.tex",
   #     showPostnikovQuiver         = false,
   #     showPlabicGraph             = false, 
   #     showPostnikovDiagram        = true, 
   #     showPostnikovDiamgramArrows = true,
   # );
end


