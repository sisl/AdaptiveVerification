using DataStructures

include("tree_utils.jl")
include("run_reluval.jl")

function adaptive_verification(root_node::NODE; min_diff = 0.1)
     q = Queue{NODE}()
     enqueue!(q, root_node)
     while length(q) > 0
        println(length(q))
        curr_node = dequeue!(q)
        next_nodes = verify_and_split!(curr_node, min_diff = min_diff)
        for i = 1:length(next_nodes)
            enqueue!(q, next_nodes[i])
        end 
     end
     return root_node
end

function verify_and_split!(curr_node::NODE; min_diff = 0.1)
    next_nodes = []
    # First run verification on possible advisories
    get_categories!(curr_node)
    # Decide if we should split
    max_diff = maximum(curr_node.ubs .- curr_node.lbs)
    split = (length(curr_node.cats) > 1) && (max_diff > min_diff)
    if split
        new_lbs, new_ubs = split_bounds(curr_node)
        for i = 1:length(new_lbs)
            push!(next_nodes, node(lbs = new_lbs[i], ubs = new_ubs[i], cats = curr_node.cats))
        end
    end
    curr_node.children = next_nodes
    return next_nodes
end

function split_bounds(curr_node::NODE; dims_to_split = [1,2])
    half_diff = (curr_node.ubs - curr_node.lbs) ./ 2
    new_lbs = Vector{Vector{Float64}}()
    new_ubs = Vector{Vector{Float64}}()

    ndims = length(dims_to_split)
    for i = 0:(2^ndims - 1)
        reg = digits(i, base = 2, pad = ndims)
        curr_lbs, curr_ubs = deepcopy(curr_node.lbs), deepcopy(curr_node.ubs)
        curr_lbs[dims_to_split] = curr_node.lbs[dims_to_split] + reg .* half_diff[dims_to_split]
        curr_ubs[dims_to_split] = curr_node.ubs[dims_to_split] - half_diff[dims_to_split] .+ reg .* half_diff[dims_to_split]
        push!(new_lbs, curr_lbs)
        push!(new_ubs, curr_ubs)
    end
    return new_lbs, new_ubs
end