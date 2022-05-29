module Postnikov

include("postnikovQuiver.jl")
export PostnikovDiagram,

    quiverFromCollection,
    constructCliqueQuiver,
    drawPostnikovDiagram


include("collections.jl")

export LabelCollection, sort!
export isNonCrossingCollection,
    isMaximalNonCrossingCollection,
    rotate!,
    isSymmetricCollection

export isNonCrossing
export whiteCliques, blackCliques

export findMaxinalNonCrossingCollections

#helper
#export isCyclicOrdered, rotateSet

include("gapqpa.jl")


include("mutation.jl")
export  mutations, mutationAtLabel


end # module
