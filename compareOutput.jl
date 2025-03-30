using JLD2
#?????????????
function arrayCompare(first, second)
    first_set = Set(eachrow(first))
    second_set = Set(eachrow(second))

    return issetequal(first_set, second_set)
end

function main()
    sizes = [1000, 3000, 10000, 30000,100000]
    setK = [2,3,5,8]
    setNumbers = 1
    resultDirectory = "result/"
    dataTypeName = ["Discs", "HollowDiscs","Square","HollowSquare","Sun","HollowSun","Circles"]
    dataType = 7
    instanceNames = Vector{String}(undef, length(sizes)*setNumbers)

    for i in 1:length(sizes)
        for j in 1:setNumbers
            k = (i-1)*setNumbers+j
            instanceNames[k] = string(dataTypeName[dataType], "_", sizes[i], "_", j)

            QHull_file = string(resultDirectory, instanceNames[k], "_QHull",".jld2")
            QHull_vertices = load(QHull_file, "V")
            
            lenk = length(setK)
            for l in 1:lenk
                alg_8ksided4CH_file = string(resultDirectory, instanceNames[k], "_8ksided4CH_", setK[l],".jld2")
                alg_8ksided4CHParallel_ver1_file = string(resultDirectory, instanceNames[k], "_8ksided4CHParallel_ver1_", setK[l],".jld2")
                alg_8ksided4CHParallel_ver2_file = string(resultDirectory, instanceNames[k], "_8ksided4CHParallel_ver2_", setK[l],".jld2")
                
                alg_8ksided4CH_vertices = load(alg_8ksided4CH_file, "V")
                alg_8ksided4CHParallel_ver1_vertices = load(alg_8ksided4CHParallel_ver1_file, "V")
                alg_8ksided4CHParallel_ver2_vertices = load(alg_8ksided4CHParallel_ver2_file, "V")
                
                println("Checking ", instanceNames[k],",k = ",setK[l])
                println("\t alg_8ksided4CH vs QHull, results are identical: ", arrayCompare(alg_8ksided4CH_vertices, QHull_vertices))
                println("\t alg_8ksided4CHParallel_ver1 vs QHull, results are identical: ", arrayCompare(alg_8ksided4CHParallel_ver1_vertices, QHull_vertices))
                println("\t alg_8ksided4CHParallel_ver2 vs QHull, results are identical: ", arrayCompare(alg_8ksided4CHParallel_ver2_vertices, QHull_vertices))
            end
        end
    end
end

main()

