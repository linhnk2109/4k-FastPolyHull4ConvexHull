using QHull

include("utils.jl")
function callQHull(points, exportCH =false, exportFile="")
        convexHull = chull(points)

    if exportCH
        exportResult(points, convexHull.vertices, exportFile)
    end

end
X= create_discs(1000)
@time callQHull(X)
