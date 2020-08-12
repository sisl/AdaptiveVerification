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
    h_inrange_pos = lbs[1] .< hsplits .< ubs[1]
    h_inrange_neg = lbs[1] .< -hsplits .< ubs[1]
    h_inrange = h_inrange_pos .| h_inrange_neg

    ḣ₀_inrange_pos = lbs[2] .< ḣ₀splits .< ubs[2]
    ḣ₀_inrange_neg = lbs[2] .< -ḣ₀splits .< ubs[2]
    ḣ₀_inrange = ḣ₀_inrange_pos .| ḣ₀_inrange_neg

    # println("hmin: $(lbs[1]), hmax: $(ubs[1])")
    # println("in range: $h_inrange")
    # println("in rangep pos: $h_inrange_pos")

    hbin = findfirst(h_inrange)
    if hbin == nothing
        hbin = findfirst(max(abs(lbs[1]), abs(ubs[1])) .< hsplits)
    end
    ḣ₀bin = findfirst(ḣ₀_inrange)
    if ḣ₀bin == nothing
        ḣ₀bin = findfirst(max(abs(lbs[2]), abs(ubs[2])) .< ḣ₀splits)
    end

    dims_to_split = []
    curr_diff = ubs .- lbs
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

# Hard-coded for 2D with splitting strategy in all dimensions
function equivalent_2D(curr_node::NODE)
    children = curr_node.children
    if length(children) == 0
        if length(curr_node.cats) == 0
            curr_node.cats = collect(0:8)
            get_categories!(curr_node, network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra03_v5_25HU_1000.nnet")
        end
        return leafnode(cats = curr_node.cats)
    elseif length(children) == 2
        # Figure out which dimension got split
        lbdiffs = children[2].lbs .- curr_node.lbs
        split_dim = findfirst(lbdiffs .!= 0.0)
        split = children[2].lbs[split_dim]

        # Create the tree
        top = kdnode(dim = split_dim, split = split)

        # Recursively add children
        top.left = equivalent_2D(children[1])
        top.right = equivalent_2D(children[2])
        return top
    else
        first_split = children[2].lbs[1]
        second_split = children[3].lbs[2]

        # Create the tree
        top = kdnode(dim = 1, split = first_split)
        top.left = kdnode(dim = 2, split = second_split)
        top.right = kdnode(dim = 2, split = second_split)

        # Recursively add children
        top.left.left = equivalent_2D(children[1])
        top.right.left = equivalent_2D(children[2])
        top.left.right = equivalent_2D(children[3])
        top.right.right = equivalent_2D(children[4])
        return top
    end
end


