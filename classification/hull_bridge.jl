## Minimal bridge used by the Section 5 Python scripts.
##
## It reads batches of 2D Float64 point sets, calls exactly
## alg_8ksided4CH or alg_8ksided4CH_parallel_ver1, and writes the returned
## vertices plus timing metadata.  Binary coordinates avoid CSV precision,
## quoting, and file-size problems at n = 10^8. No non-standard Julia package
## is required.
##
## ARGS:
##   1. points.f64   raw little-endian Float64 values: x1,y1,x2,y2,...
##   2. manifest.tsv id, group, 1-based point range, k, mode, timing settings
##   3. hulls.f64    raw returned vertices: x1,y1,x2,y2,...
##   4. results.tsv  vertex ranges, timing runs, and counts

using Statistics

include("8ksidedPolygon4CH.jl")
include("8ksidedPolygon4CH_parallel_ver1.jl")


struct Job
    id::String
    group::String
    start_row::Int
    end_row::Int
    k::Int
    mode::String
    time_it::Bool
    repeats::Int
end


function read_points_binary(path)
    bytes = read(path)
    length(bytes) % 16 == 0 || error("points file size must be a multiple of 16 bytes")
    values = reinterpret(Float64, bytes)
    return permutedims(reshape(values, 2, :))
end


function read_manifest(path)
    lines = readlines(path)
    isempty(lines) && error("empty manifest")
    expected_header = "id\tgroup\tstart_row\tend_row\tk\tmode\ttime_it\trepeats"
    lines[1] == expected_header || error("unexpected manifest header: $(lines[1])")
    jobs = Job[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, '\t'; keepempty=true)
        length(fields) == 8 || error("manifest row must contain 8 tab-separated fields: $line")
        job = Job(
            fields[1],
            fields[2],
            parse(Int, fields[3]),
            parse(Int, fields[4]),
            parse(Int, fields[5]),
            fields[6],
            fields[7] == "1",
            parse(Int, fields[8]),
        )
        job.k >= 1 || error("job $(job.id): k must be >= 1")
        job.mode in ("seq", "par") || error("job $(job.id): invalid mode $(job.mode)")
        job.repeats >= 1 || error("job $(job.id): repeats must be >= 1")
        push!(jobs, job)
    end
    return jobs
end


function run_hull(points, k, mode)
    n = size(points, 1)
    if n == 0
        return Vector{Vector{Float64}}()
    elseif n == 1
        return [[Float64(points[1, 1]), Float64(points[1, 2])]]
    elseif n == 2
        first_point = [Float64(points[1, 1]), Float64(points[1, 2])]
        second_point = [Float64(points[2, 1]), Float64(points[2, 2])]
        return first_point == second_point ? [first_point] : [first_point, second_point]
    end

    raw = if mode == "seq"
        alg_8ksided4CH(points, k)
    elseif mode == "par"
        alg_8ksided4CH_parallel_ver1(points, k)
    else
        error("unknown mode: $mode")
    end
    return [[Float64(point[1]), Float64(point[2])] for point in raw]
end


function validate_ranges(jobs, number_of_points)
    ids = Set{String}()
    for job in jobs
        job.id in ids && error("duplicate job id: $(job.id)")
        push!(ids, job.id)
        number_for_job = max(0, job.end_row - job.start_row + 1)
        if number_for_job == 0
            job.start_row == job.end_row + 1 || error("job $(job.id): invalid empty range")
            1 <= job.start_row <= number_of_points + 1 || error("job $(job.id): range is out of bounds")
        else
            1 <= job.start_row <= job.end_row <= number_of_points ||
                error("job $(job.id): range is out of bounds")
        end
    end
end


function median_seconds(values)
    return isempty(values) ? nothing : median(values)
end


function main()
    length(ARGS) == 4 || error(
        "usage: julia --threads 2 hull_bridge.jl points.f64 manifest.tsv hulls.f64 results.tsv"
    )
    points_path, manifest_path, hulls_path, results_path = ARGS

    all_points = read_points_binary(points_path)
    jobs = read_manifest(manifest_path)
    validate_ranges(jobs, size(all_points, 1))
    println(
        "hull_bridge.jl: ", size(all_points, 1), " stored points, ", length(jobs),
        " jobs, Threads.nthreads()=", Threads.nthreads(),
    )

    # Compile both code paths before any reported timing.
    warmup = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0; 0.5 0.5; 2.0 0.3]
    run_hull(warmup, 1, "seq")
    run_hull(warmup, 1, "par")
    # Real jobs are SubArray views into the shared binary matrix.  Julia
    # specializes on argument type, so warm that code path too; otherwise the
    # first reported timing run includes a second compilation.
    warmup_view = @view warmup[:, :]
    run_hull(warmup_view, 1, "seq")
    run_hull(warmup_view, 1, "par")

    inputs = Dict{String, Any}()
    n_inputs = Dict{String, Int}()
    for job in jobs
        if job.end_row < job.start_row
            inputs[job.id] = @view all_points[1:0, :]
            n_inputs[job.id] = 0
        else
            inputs[job.id] = @view all_points[job.start_row:job.end_row, :]
            n_inputs[job.id] = job.end_row - job.start_row + 1
        end
    end

    hulls = Dict{String, Any}()
    timing_runs = Dict{String, Vector{Float64}}()
    timed_groups = Dict{String, Vector{Job}}()
    group_order = String[]

    for job in jobs
        if job.time_it
            if !haskey(timed_groups, job.group)
                timed_groups[job.group] = Job[]
                push!(group_order, job.group)
            end
            push!(timed_groups[job.group], job)
        else
            hulls[job.id] = run_hull(inputs[job.id], job.k, job.mode)
            timing_runs[job.id] = Float64[]
        end
    end

    for (group_index, group_name) in enumerate(group_order)
        group_jobs = timed_groups[group_name]
        repeats = group_jobs[1].repeats
        all(job -> job.repeats == repeats, group_jobs) ||
            error("timed group $group_name mixes repeat counts")
        # The real view has a UnitRange row index (different from the full
        # warmup view's Slice type), so compile one representative call before
        # starting the clock.  This is an unreported JIT warmup, not a repeat.
        warm_group_hulls = Dict{String, Any}()
        for job in group_jobs
            warm_group_hulls[job.id] = run_hull(inputs[job.id], job.k, job.mode)
        end
        runs = Float64[]
        final_hulls = Dict{String, Any}()
        for _ in 1:repeats
            start_ns = time_ns()
            for job in group_jobs
                final_hulls[job.id] = run_hull(inputs[job.id], job.k, job.mode)
            end
            push!(runs, (time_ns() - start_ns) / 1.0e9)
        end
        for job in group_jobs
            hulls[job.id] = final_hulls[job.id]
            timing_runs[job.id] = copy(runs)
        end
        println(
            "hull_bridge.jl: timed group ", group_index, "/", length(group_order),
            " (", group_name, ") median=", round(median(runs); digits=6), "s",
        )
    end

    vertex_ranges = Dict{String, Tuple{Int, Int}}()
    next_vertex = 1
    open(hulls_path, "w") do output
        for job in jobs
            hull = hulls[job.id]
            start_vertex = next_vertex
            for point in hull
                write(output, Float64(point[1]))
                write(output, Float64(point[2]))
                next_vertex += 1
            end
            vertex_ranges[job.id] = (start_vertex, next_vertex - 1)
        end
    end

    open(results_path, "w") do output
        println(output, "id\ttime_s\ttime_runs_s\tn_verts\tn_input\tstart_vertex\tend_vertex")
        for job in jobs
            runs = timing_runs[job.id]
            median_value = median_seconds(runs)
            median_text = median_value === nothing ? "" : string(median_value)
            runs_text = join(runs, ";")
            start_vertex, end_vertex = vertex_ranges[job.id]
            println(
                output,
                job.id, '\t', median_text, '\t', runs_text, '\t', length(hulls[job.id]), '\t',
                n_inputs[job.id], '\t', start_vertex, '\t', end_vertex,
            )
        end
    end
    println("hull_bridge.jl: wrote ", hulls_path, " and ", results_path)
end


main()
