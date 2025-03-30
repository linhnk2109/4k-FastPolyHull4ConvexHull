using LinearAlgebra
using Random
using LoopVectorization
include("utils.jl")

@inline function starting_vertices(X,k)
    Is = [[] for _ in 1:(4k+1)]
    len = size(X,1)
    for j in 1:2k
        d = pi*(j-1)/2k
        e_x = cos(d)
        e_y = sin(d)
        max_j = e_x*X[1,1]+e_y*X[1,2]
        min_j = e_x*X[1,1]+e_y*X[1,2]
        if j >0 && j <=k
            p_max = [X[1,1],X[1,2]]
            p_min = [X[1,1],X[1,2]]
            
            p_max2 = [X[1,1],X[1,2]]
            p_min2 = [X[1,1],X[1,2]]

            @inbounds for i in 2:len
                x = X[i,1]
                y = X[i,2]
                t = e_x*x+ e_y*y
                #I^j
                max_j_eq_t = max_j == t
                max_j_lt_t = max_j < t

                y_gt_y_max = y > p_max[2]
                y_lt_y_min = y < p_min[2]

                max_j = max_j_lt_t ? t : max_j
                p_max[1] = max_j_lt_t ? x : p_max[1]
                p_min[1] = max_j_lt_t ? x : p_min[1]
                p_max[2] = max_j_lt_t ? y : p_max[2]
                p_min[2] = max_j_lt_t ? y : p_min[2]

                p_min[1] = (max_j_eq_t & y_lt_y_min) ? x : p_min[1]
                p_min[2] = (max_j_eq_t & y_lt_y_min) ? y : p_min[2]
                p_max[1] = (max_j_eq_t & y_gt_y_max) ? x : p_max[1]
                p_max[2] = (max_j_eq_t & y_gt_y_max) ? y : p_max[2]

                #I^(j+2k)
                min_j_eq_t = min_j == t 
                min_j_gt_t =  min_j > t
                
                y_gt_y_max2 = y > p_max2[2]
                y_lt_y_min2 = y < p_min2[2]

                min_j = min_j_gt_t ? t : min_j
                p_min2[1] = min_j_gt_t ? x : p_min2[1]
                p_max2[1] = min_j_gt_t ? x : p_max2[1]
                p_min2[2] = min_j_gt_t ? y : p_min2[2]
                p_max2[2] = min_j_gt_t ? y : p_max2[2]

                p_min2[1] = (min_j_eq_t & y_lt_y_min2) ? x : p_min2[1]
                p_min2[2] = (min_j_eq_t & y_lt_y_min2) ? y : p_min2[2]
                p_max2[1] = (min_j_eq_t & y_gt_y_max2) ? x : p_max2[1]
                p_max2[2] = (min_j_eq_t & y_gt_y_max2) ? y : p_max2[2]
        
            end
            Is[j] = [p_min,p_max]
            Is[j+2k] = [p_max2,p_min2]
        else
            p_max = [X[1,1],X[1,2]]
            p_min = [X[1,1],X[1,2]]

            p_max2 = [X[1,1],X[1,2]]
            p_min2 = [X[1,1],X[1,2]]
            @inbounds for i in 2:len
                x = X[i,1]
                y = X[i,2]
                t = e_x*x+ e_y*y
                #I^j
                max_j_eq_t = max_j == t
                max_j_lt_t = max_j < t

                x_gt_x_max = x > p_max[1]
                x_lt_x_min = x < p_min[1]

                max_j = max_j_lt_t ? t : max_j
                p_max[1] = max_j_lt_t ? x : p_max[1]
                p_min[1] = max_j_lt_t ? x : p_min[1]
                p_max[2] = max_j_lt_t ? y : p_max[2]
                p_min[2] = max_j_lt_t ? y : p_min[2]

                p_min[1] = (max_j_eq_t & x_lt_x_min) ? x : p_min[1]
                p_min[2] = (max_j_eq_t & x_lt_x_min) ? y : p_min[2]
                p_max[1] = (max_j_eq_t & x_gt_x_max) ? x : p_max[1]
                p_max[2] = (max_j_eq_t & x_gt_x_max) ? y : p_max[2]

                #I^(j+2k)
                min_j_eq_t = min_j == t 
                min_j_gt_t =  min_j > t
                
                x_gt_x_max2 = x > p_max2[1]
                x_lt_x_min2 = x < p_min2[1]

                min_j = min_j_gt_t ? t : min_j
                p_min2[1] = min_j_gt_t ? x : p_min2[1]
                p_max2[1] = min_j_gt_t ? x : p_max2[1]
                p_min2[2] = min_j_gt_t ? y : p_min2[2]
                p_max2[2] = min_j_gt_t ? y : p_max2[2]

                p_min2[1] = (min_j_eq_t & x_lt_x_min2) ? x : p_min2[1]
                p_min2[2] = (min_j_eq_t & x_lt_x_min2) ? y : p_min2[2]
                p_max2[1] = (min_j_eq_t & x_gt_x_max2) ? x : p_max2[1]
                p_max2[2] = (min_j_eq_t & x_gt_x_max2) ? y : p_max2[2]
            end
            Is[j] = [p_max,p_min]
            Is[j+2k] = [p_min2,p_max2]
        end        
    end
    Is[4k+1] = [Is[1][1]]
    return Is
end

@inline function find_P(X, I,k)
    Pp = [zeros(0,2) for _ in 1:4k]
    len = size(X,1)
    epsilon = 10^-8
    for j in 1:4k
        alpha = I[j+1][1][2] - I[j][2][2]
        beta = I[j+1][1][1] - I[j][2][1]
        gamma = alpha*I[j][2][1] - beta*I[j][2][2]
        mask = Vector{Bool}(undef,len)
        @inbounds @avx for i in 1:len
            mask[i] = X[i,1]*alpha - X[i,2]*beta > gamma 
        end
        Pp[j] = vcat(X[mask,:],[I[j][2][1] I[j][2][2]; I[j+1][1][1] I[j+1][1][2]])
    end
    return Pp
end

@inline function reduce_1(P)
    P = sortslices(P, dims = 1, by = x -> -x[1])
    Pp = [P[1,:]]
    x_min = P[1,1]
    y_max = P[1,2]
    len = size(P,1)
    @inbounds for i in 2:len
        x = P[i,1]
        y = P[i,2]
        if y > y_max
            if x != x_min
                push!(Pp,P[i,:])
                x_min = x
                y_max = y
            else
                Pp[end] = P[i,:]
                y_max = y
            end
        end
    end
    return Pp
end

@inline function reduce_2(P)
    P = sortslices(P, dims = 1, by = x -> x[1])
    Pp = [P[1,:]]
    x_max = P[1,1]
    y_max = P[1,2]
    len = size(P,1)
    @inbounds for i in 2:len
        x = P[i,1]
        y = P[i,2]
        if y > y_max
            if x != x_max
                push!(Pp,P[i,:])
                x_max = x
                y_max = y
            else
                Pp[end] = P[i,:]
                y_max = y
            end
        end
    end
    return Pp
end

@inline function reduce_3(P)
    P = sortslices(P, dims = 1, by = x -> x[1])
    Pp = [P[1,:]]
    x_max = P[1,1]
    y_min = P[1,2]
    len = size(P,1)
    @inbounds for i in 2:len
        x = P[i,1]
        y = P[i,2]
        if y < y_min 
            if x != x_max
                push!(Pp,P[i,:])
                x_max = x
                y_min = y
            else
                Pp[end] = P[i,:]
                y_min = y
            end
        end
    end
    return Pp
end

@inline function reduce_4(P)
    P = sortslices(P, dims = 1, by = x -> -x[1])
    Pp = [P[1,:]]
    x_min = P[1,1]
    y_min = P[1,2]
    len = size(P,1)
    @inbounds for i in 2:len
        x = P[i,1]
        y = P[i,2]
        if y < y_min 
            if x != x_min
                push!(Pp,P[i,:])
                x_min = x
                y_min = y
            else
                Pp[end] = P[i,:]
                y_min = y
            end 
        end
    end
    return Pp
end

@inline function algorithm_1(P)
    l = length(P)
    l<3 && return P
    V = Vector{Vector{Float64}}(undef, ceil(Int, l/10))
    sizeV = 1
    V[1] = P[1]

    stack_size = 1
    i =2
    
    @inbounds while i <= l
        while stack_size >=2 && orient(V[stack_size-1],V[stack_size],P[i])<=0 stack_size-=1 end
        sizeV = stack_size
        sizeV+=1
        if sizeV > length(V)
            resize!(V, ceil(Int, length(V) * 1.5))
        end
        V[sizeV] = P[i]
        stack_size = sizeV
        i += 1
    end
    resize!(V, sizeV)
    return V
end

@inline function algorithm_2(P)
    l = length(P)
    l<3 && return P
    V = Vector{Vector{Float64}}(undef, ceil(Int, l/10))
    sizeV = 1
    V[1] = P[1]

    stack_size = 1
    i =2
    
    @inbounds while i <= l
        while stack_size >=2 && orient(V[stack_size-1],V[stack_size],P[i])>=0 stack_size-=1 end
        sizeV = stack_size
        sizeV+=1
        if sizeV > length(V)
            resize!(V, ceil(Int, length(V) * 1.5))
        end
        V[sizeV] = P[i]
        stack_size = sizeV
        i += 1
    end
    resize!(V, sizeV)
    return V
end

@inline function algorithm_3(P)
    l = length(P)
    l<3 && return P
    V = Vector{Vector{Float64}}(undef, ceil(Int, l/10))
    sizeV = 1
    V[1] = P[1]

    stack_size = 1
    i =2
    
    @inbounds while i <= l
        while stack_size >=2 && orient(V[stack_size-1],V[stack_size],P[i])<=0 stack_size-=1 end
        sizeV = stack_size
        sizeV+=1
        if sizeV > length(V)
            resize!(V, ceil(Int, length(V) * 1.5))
        end
        V[sizeV] = P[i]
        stack_size = sizeV
        i += 1
    end
    resize!(V, sizeV)
    return V
end

@inline function algorithm_4(P)
    l = length(P)
    l<3 && return P
    V = Vector{Vector{Float64}}(undef, ceil(Int, l/10))
    sizeV = 1
    V[1] = P[1]

    stack_size = 1
    i =2
    
    @inbounds while i <= l
        while stack_size >=2 && orient(V[stack_size-1],V[stack_size],P[i])>=0 stack_size-=1 end
        sizeV = stack_size
        sizeV+=1
        if sizeV > length(V)
            resize!(V, ceil(Int, length(V) * 1.5))
        end
        V[sizeV] = P[i]
        stack_size = sizeV
        i += 1
    end
    resize!(V, sizeV)
    return V
end
@inline function alg_8ksided4CH(X,k)
    Is = starting_vertices(X,k)
    P = find_P(X,Is,k)
    V = [[] for _ in 1:4k]

    for i in 1:4k
        if i >0 && i <=k
            Pp = reduce_1(P[i])
            V[i] = algorithm_1(Pp)
        elseif i >k && i<=2k
            Pp = reduce_2(P[i])
            V[i] = reverse(algorithm_2(Pp))
        elseif i >2k && i <=3k
            Pp = reduce_3(P[i])
            V[i] = algorithm_3(Pp)
        else
            Pp = reduce_4(P[i])
            V[i] = reverse(algorithm_4(Pp))
        end   
    end
    hull = V[1]
    for i in 2:4k 
        if hull[end] == V[i][1]
            append!(hull, V[i][2:end])
        else
            append!(hull, V[i])
        end
    end
    if hull[1] == hull[end]
        pop!(hull)
    end
    return hull
end


@inline function alg_8ksided4CH(X,k, exportFile)
    hull = alg_8ksided4CH(X,k)
    exportResult(hull, exportFile)
end

X = create_discs(1000)
hull = @time alg_8ksided4CH(X,2)
