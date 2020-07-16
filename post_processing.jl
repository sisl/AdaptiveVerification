hsplits = [150, 800, 1250, 4000, Inf]
hwidths = [25, 50, 100, 250, 1000]

ḣ₀splits = [30, 50, Inf]
ḣ₀widths = [3, 5, 10]

# ḣ₁splits = [0, 30, 50]
# ḣ₁widths = [3, 5, 10]

function discretize_for_mdp(root_node::NODE)
    q = Queue{NODE}()
    enqueue!(q, root_node)
    while length(q) > 0
        curr_node = dequeue!(q)
        next_nodes = length(curr_node.children) > 0 ? curr_node.children : discretize!(curr_node)
        for i = 1:length(next_nodes)
            enqueue!(q, next_nodes[i])
        end 
    end
end

function discretize!(curr_node::NODE)
    ubs = unnormalize_bounds(curr_node.ubs)
    lbs = unnormalize_bounds(curr_node.lbs)
    
    # Decide if need to split
    curr_diff = ubs .- lbs
    hbin = findfirst(max(abs(lbs[1]), abs(ubs[1])) .< hsplits)
    ḣ₀bin = findfirst(max(abs(lbs[2]), abs(ubs[2])) .< ḣ₀splits)
    dims_to_split = []
    curr_diff[1] > hwidths[hbin] ? push!(dims_to_split, 1) : nothing
    curr_diff[2] > ḣ₀widths[ḣ₀bin] ? push!(dims_to_split, 2) : nothing

    next_nodes = []

    if length(dims_to_split) > 0
        # println(curr_diff)
        new_lbs, new_ubs = split_bounds(curr_node, dims_to_split = dims_to_split)
        for i = 1:length(new_lbs)
            push!(next_nodes, node(lbs = new_lbs[i], ubs = new_ubs[i], cats = curr_node.cats))
        end
        curr_node.children = next_nodes
    end

    return next_nodes
end

"""
KD Trees
"""

# mutable struct KDNODE
#     dim::Int32
#     split::Float64
#     left::KDNODE
#     right::KDNODE
# end

function equivalent_KD(curr_node::NODE)
    children = curr_node.children
    lbs = curr_node.lbs
    ubs = curr_node.ubs

    lb_diffs = [child.lbs .- lbs for child in children]

    # Find out where the splits are
    ndims = length(lbs)
    dims = []
    splits = []
    for i = 1:length(lb_diffs)
        diff_dims = findall(lb_diffs[i] .> 0.0)
        for dim in diff_dims
            if !(dim in dims)
                push!(dims, dim)
                push!(splits, children[i].lbs[dim])
            end
        end
    end
end

