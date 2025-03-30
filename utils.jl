using Random
using Printf
using CSV
using DataFrames
using StatsBase
using JLD2
function create_discs(n)
    points = Matrix{Float64}(undef,n,2)
    alpha = rand(0:359,n)
    bankinh = 100*sqrt.(rand(n))
    points[:,1] = bankinh .* sin.(alpha)
    points[:,2] = bankinh .* cos.(alpha)
    # points = round.(points,digits=6)
    return points
end

function create_hollowDiscs(n)
    points = Matrix{Float64}(undef,n,2)
    alpha = rand(0:359,n)
    bankinh = 75 .+ 25 .* sqrt.(rand(n))
    points[:,1] = bankinh .* sin.(alpha)
    points[:,2] = bankinh .* cos.(alpha)
    # points = round.(points,digits=1)
    return points
end

function create_square(n)
    points = Matrix{Float64}(undef,n,2)
    points[:,1] = rand(-100:100,n)
    points[:,2] = rand(-100:100,n)
    # points = round.(points,digits=1)
    return points
end

function create_hollowSquare(n)
    points = Matrix{Float64}(undef,n,2)
    for i in 1:n 
        if i <= n÷4
            points[i,1] = 75 + 25 * rand()
            points[i,2] = 100 * rand()
        elseif i <= n÷2
            points[i,1] = 25 * rand()
            points[i,2] = 100 * rand()
        elseif i <= 3n÷4
            points[i,1] = 100 * rand()
            points[i,2] = 75 + 25 * rand()
        else
            points[i,1] = 100 * rand()
            points[i,2] = 25 * rand()
        end
    end
    # points = round.(points,digits=1)
    return points
end

function create_sun(n)
    points = Matrix{Float64}(undef,n,2)
    alpha = rand(0:99,n)
    bankinh = rand(1:300,n)
    points[:,1] = bankinh .* sin.(alpha)
    points[:,2] = bankinh .* cos.(alpha)
    # points = round.(points,digits=1)
    return points
end

function create_hollowSun(n)
    points = Matrix{Float64}(undef,n,2)
    alpha = rand(0:99,n)
    bankinh = 380 .+ rand(1:300,n)
    points[:,1] = bankinh .* sin.(alpha)
    points[:,2] = bankinh .* cos.(alpha)
    # points = round.(points,digits=1)
    return points
end

function create_circles(n)
    points = Matrix{Float64}(undef,n,2)
    scale1 = rand(n÷2)
    scale2 = rand(n-n÷2)
    points[1:n÷2,1] = 200 .* scale1 .-100
    points[1:n÷2,2] = sqrt.(10000 .- points[1:n÷2,1].^2)
    points[n÷2+1:n,1] = 200 .* scale2 .-100
    points[n÷2+1:n,2] = -sqrt.(10000 .- points[n÷2+1:n,1].^2)
    # points = round.(points,digits=6)
    return points
end
@inline function orient(p1, p2, p3)
    dx1 = p2[1] - p1[1]
    dy1 = p2[2] - p1[2]
    dx2 = p3[1] - p1[1]
    dy2 = p3[2] - p1[2]
    return dx1 * dy2 - dy1 * dx2 
end
function report(bm, ignore::Int)
    # times in nano seconds hence multiplied by 1.e-9
    sorted_times = sort(bm.times)
    s = max(1, length(sorted_times)-ignore)
    min_time = minimum(sorted_times[1:s])*1.e-9
    max_time = maximum(sorted_times[1:s])*1.e-9
    mean_time = mean(sorted_times[1:s])*1.e-9
    geomean_time = geomean(sorted_times[1:s])*1.e-9
    median_time = median(sorted_times[1:s])*1.e-9
    runs = length(bm)
    @printf("running time (seconds):  min=%.5f", min_time)
    @printf(" max=%.5f", max_time)
    @printf(" mean=%.5f", mean_time)
    @printf(" geomean=%.5f", geomean_time)
    @printf(" median=%.5f", median_time)
    print(" runs=", runs)
    println(" ignore=", ignore)
    return [min_time, max_time, mean_time, geomean_time, median_time, runs, ignore]
end

function exportReport(names, runningTime, file_name)
    df = DataFrame(instance = names,
                    min = map(i -> runningTime[i,1], collect(1:length(names))),
                    max = map(i -> runningTime[i,2], collect(1:length(names))),
                    mean = map(i -> runningTime[i,3], collect(1:length(names))),
                    geomean = map(i -> runningTime[i,4], collect(1:length(names))),
                    median = map(i -> runningTime[i,5], collect(1:length(names))),
                    runs = map(i -> runningTime[i,6], collect(1:length(names))),
                    ignore = map(i -> runningTime[i,7], collect(1:length(names))))
    CSV.write(file_name, df)
end

function exportResult(points, vertices, exportFile)
    V = Vector{Vector{Float64}}(undef,length(vertices))
    for i in 1:length(vertices) 
        V[i] = points[vertices[i],:]
    end
    jldsave(string(exportFile, ".jld2"); V)
    # df = DataFrame(index = vertices,
    #                 x = map(v -> points[v,1], vertices),
    #                 y = map(v -> points[v,2], vertices))
    # CSV.write(string(exportFile, ".csv"), df)
end

function exportResult(points, exportFile)
    V = Vector{Vector{Float64}}(undef,length(points))
    for i in 1:length(points)
        V[i] = points[i]
    end
    jldsave(string(exportFile, ".jld2"); V)
    # df = DataFrame(V)
    # CSV.write(string(exportFile, ".csv"), df)
end

@inline function mySplit(startIndex, endIndex, numberOfSet)
    step = Int(ceil((endIndex-startIndex+1)/numberOfSet))
    grid = [startIndex:step:endIndex;endIndex+1]
    return grid
end