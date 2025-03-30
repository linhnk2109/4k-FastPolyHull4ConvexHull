using Random
using BenchmarkTools
using Distributions
using CGAL
using IterTools
using TimerOutputs
include("utils.jl")
include("8ksidedPolygon4CH.jl")
include("8ksidedPolygon4CH_parallel_ver1.jl")
#include("8ksidedPolygon4CH_parallel_ver2.jl")
#include("quickHull.jl")

function main()
    #compareOutput
    # sizes = [1000, 3000, 10000, 30000,100000]
    # sizes = [1000, 3000, 10000, 30000,100000, 300000, 1000000] #test output
    sizes = [10_000, 31_600, 100_000, 316_000, 1_000_000, 3_160_000, 10_000_000, 31_600_000, 100_000_000] 
    setK = [6]
    setNumbers = 1
    resultDirectory = "result/"
    # set benchmarking = true if want to benchmark
    benchmarking = true
    # benchmarking = false
    # only export results if exportResult = true and benchmarking = false
    exportResult = true
    # dataType: 
        #1 for discs type, 
        #2 for hollow discs type
        #3 for square type
        #4 for hollow square type
        #5 for sun type
        #6 for hollow sun type
        #7 for circles type
    dataType = 3
    dataTypeName = ["Discs", "HollowDiscs","Square","HollowSquare","Sun","HollowSun","Circles"]
    
    instanceNames = Vector{String}(undef, length(sizes)*setNumbers)

    runningTimeQHull = Matrix{Float64}(undef, length(instanceNames),7) 
    runningTime8ksided4CH = [Matrix{Float64}(undef, length(instanceNames),7) for _ in 1:length(setK)]
    runningTime8ksided4CHParallel_ver1 = [Matrix{Float64}(undef, length(instanceNames),7) for _ in 1:length(setK)]
    runningTime8ksided4CHParallel_ver2 = [Matrix{Float64}(undef, length(instanceNames),7) for _ in 1:length(setK)]
    runningTimeAklToussaint = Matrix{Float64}(undef, length(instanceNames), 7)
    runningTimeBykat = Matrix{Float64}(undef, length(instanceNames), 7)
    runningTimeEddy = Matrix{Float64}(undef, length(instanceNames), 7)
    runningTimeGrahamAndrew = Matrix{Float64}(undef, length(instanceNames), 7)
    runningTimeJavisMarch = Matrix{Float64}(undef, length(instanceNames), 7)
    runningTimeConvexHull2 = Matrix{Float64}(undef, length(instanceNames), 7)
    runs = 5
    ignore = 2
    Random.seed!(42)
    for i in 1:length(sizes)
        for j in 1:setNumbers
            k = (i-1)*setNumbers+j
            instanceNames[k] = string(dataTypeName[dataType], "_", sizes[i], "_", j)

            println()
            println("Consider instance ", instanceNames[k])

            # create random data
            points = Matrix{Float64}(undef, sizes[i], 2)
            if dataType == 1
                points = create_discs(sizes[i])
            elseif dataType == 2
                points = create_hollowDiscs(sizes[i])
            elseif dataType == 3
                points = create_square(sizes[i])
            elseif dataType == 4
                points = create_hollowSquare(sizes[i])
            elseif dataType == 5
                points = create_sun(sizes[i])
            elseif dataType == 6
                points = create_hollowSun(sizes[i])
            else
                points = create_circles(sizes[i])
            end
            points_CGAL = Vector{Point2}(undef, size(points, 1))

            for t in 1:length(points_CGAL)
                points_CGAL[t] = Point2(points[t,1],points[t,2])
            end
            if benchmarking
                
                lenk = length(setK)
                for l in 1:lenk
                    K_ = setK[l]
                    #=
                    println("8k_Sided_4_OCH; k = ",K_)
                    bm8ksided4CH = run(@benchmarkable alg_8ksided4CH($points,$K_) samples=runs seconds=10000)
                    runningTime8ksided4CH[l][k,:] = report(bm8ksided4CH, ignore)
=#
                    println("8k_Sided_4_CH_Parallel_ver1; k = ",K_)
                    bm8ksided4CHParallel_ver1 = run(@benchmarkable alg_8ksided4CH_parallel_ver1($points,$K_) samples=runs seconds=10000)
                    runningTime8ksided4CHParallel_ver1[l][k,:] = report(bm8ksided4CHParallel_ver1, ignore)
#=
                    println("8k_Sided_4_CH_Parallel_ver2; k = ",K_)
                    bm8ksided4CHParallel_ver2 = run(@benchmarkable alg_8ksided4CH_parallel_ver2($points,$K_) samples=runs seconds=10000)
                    runningTime8ksided4CHParallel_ver2[l][k,:] = report(bm8ksided4CHParallel_ver2, ignore)
                    =#
                end
                #=
                println("QHull")
                bmQHull = run(@benchmarkable callQHull($points) samples=runs seconds=10000)
                runningTimeQHull[k,:] = report(bmQHull, ignore)
                
                println("Javis March")
                bmQH = run(@benchmarkable ch_jarvis($points_CGAL) samples= runs seconds=10000)
                runningTimeJavisMarch[k,:] = report(bmQH, ignore)
                
                println("Convex Hull 2")
                bmQH = run(@benchmarkable convex_hull_2($points_CGAL) samples= runs seconds=10000)
                runningTimeConvexHull2[k,:] = report(bmQH, ignore)
                
                println("Akl Toussaint heuristic")
                bmAklToussaint = run(@benchmarkable ch_akl_toussaint($points_CGAL) samples=runs seconds=10000)
                runningTimeAklToussaint[k,:] = report(bmAklToussaint, ignore)

                println("Bykat-quickHull-non-recursive version")
                bmBykat = run(@benchmarkable ch_bykat($points_CGAL) samples=runs seconds=10000)
                runningTimeBykat[k,:] = report(bmBykat, ignore)
                
                println("Eddy - a version of quickHull algorithm")
                bmEddy = run(@benchmarkable ch_eddy($points_CGAL) samples=runs seconds=10000)
                runningTimeEddy[k,:] = report(bmEddy, ignore)

                println("Graham Andrew - Graham' scan")
                bmGrahamAndrew = run(@benchmarkable ch_graham_andrew($points_CGAL) samples=runs seconds=10000)
                runningTimeGrahamAndrew[k,:] = report(bmGrahamAndrew, ignore)
                =#
            else
                lenk = length(setK)
                for l in 1:lenk
                    K_ = setK[l]
                    exportFile8ksided4CH = string(resultDirectory, instanceNames[k], "_8ksided4CH_", setK[l])
                    alg_8ksided4CH(points,K_,exportFile8ksided4CH)

                    exportFile8ksided4CHParallel_ver1 = string(resultDirectory, instanceNames[k], "_8ksided4CHParallel_ver1_", setK[l])
                    alg_8ksided4CH_parallel_ver1(points,K_,exportFile8ksided4CHParallel_ver1)

                    exportFile8ksided4CHParallel_ver2 = string(resultDirectory, instanceNames[k], "_8ksided4CHParallel_ver2_", setK[l])
                    alg_8ksided4CH_parallel_ver2(points,K_,exportFile8ksided4CHParallel_ver2)
                end
                
                exportFileQHull = string(resultDirectory, instanceNames[k], "_QHull")
                callQHull(points, exportResult,exportFileQHull)
            end
        end
    end

    if benchmarking
        baseName = string(resultDirectory, dataTypeName[dataType], "_")
        
        for l in 1:length(setK)
            #exportReport(instanceNames, runningTime8ksided4CH[l], string(baseName, "8ksided4CH_", setK[l], "_running_time.csv"))
            exportReport(instanceNames, runningTime8ksided4CHParallel_ver1[l], string(baseName, "8ksided4CHParallel_ver1_", setK[l], "_running_time.csv"))
            #exportReport(instanceNames, runningTime8ksided4CHParallel_ver2[l], string(baseName, "8ksided4CHParallel_ver2_", setK[l], "_running_time.csv"))
            
        end
#=
        exportReport(instanceNames, runningTimeQHull, string(baseName,"QHull_running_time.csv"))
        
        exportReport(instanceNames, runningTimeJavisMarch, string(baseName,"JavisMarch_running_time.csv"))
        exportReport(instanceNames, runningTimeConvexHull2, string(baseName,"ConvexHull2_running_time.csv"))
        exportReport(instanceNames, runningTimeAklToussaint, string(baseName,"Akl_Toussaint_running_time.csv"))
        exportReport(instanceNames, runningTimeBykat, string(baseName,"Bykat_running_time.csv"))
        exportReport(instanceNames, runningTimeEddy, string(baseName,"Eddy_running_time.csv"))
        exportReport(instanceNames, runningTimeGrahamAndrew, string(baseName,"Graham_Andrew_running_time.csv"))
 =#       
    end
end

main()