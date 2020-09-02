hsplits = [150, 800, 1250, 4000, Inf]
hwidths = [25, 50, 100, 250, 1000]

ḣ₀splits = [30, 50, Inf]
ḣ₀widths = [3, 5, 10]

ḣ₁splits = [30, 50, Inf]
ḣ₁widths = [3, 5, 10]

"""
Interior nodes:
nodes (Int64): left index (right index is left index + 1 always!)
dims (BitArray): one hot encoding of dimension
splits (Float64): split values

Leaf nodes:
nodes (Int64): -leaf_data index (know it is leaf node if this is a negative value)
dims (BitArray): N/A will just keep as garbage value or zero
splits (Float64): N/A will just keep as garbage value or zero
leaf_data (BitArray): one hot encoding of possible advisories
"""

function write_to_arrays(tree::VERIFYNODE; max_dim = 2)
    # First, get the number of interior and leaf nodes
    ninterior = num_interior(tree)
    nleaves = num_leaves(tree)
    ntot = ninterior + nleaves

    # Initialize arrays to store node integers and floats and leaf data
    nodes = zeros(Int32, ntot)
    splits = zeros(Float64, ntot)
    # Initializes bit array to all false (going to use one hot encoding)
    dims = falses(max_dim, ntot)
    leaf_data = falses(9, nleaves)

    # Start everything at first index
    node_index = 1
    leaf_index = 1

    # Start traversing the tree and filling in the arrays
    s = Stack{VERIFYNODE}()
    ind_s = Stack{Int32}()
    push!(s, tree)
    push!(ind_s, node_index)

    while length(s) > 0
        curr = pop!(s)
        curr_index = pop!(ind_s)

        if typeof(curr.left) == LEAFNODE # is leaf
            # Fill in everything for leaves
            nodes[curr_index] = -leaf_index
            inds = curr.cats .+ 1 # cats start at zero
            leaf_data[inds, leaf_index] .= true
            leaf_index += 1
        else # is interior
            # Fill in everything for interior
            dims[curr.dim, curr_index] = true
            splits[curr_index] = curr.split

            nodes[curr_index] = node_index + 1

            # Move to the next nodes
            push!(s, curr.left)
            push!(s, curr.right)
            push!(ind_s, node_index + 1)
            push!(ind_s, node_index + 2)
            node_index += 2
        end
    end
    return TREE_ARRAYS(nodes, splits, dims, leaf_data)
end

function write_to_files(ta::TREE_ARRAYS, prefix)
    # Nodes file
    n = open("$(prefix)_n.bin", "w+")
    write(n, length(ta.nodes))
    write(n, ta.nodes)
    close(n)

    # Splits file
    s = open("$(prefix)_s.bin", "w+")
    write(s, length(ta.splits))
    write(s, ta.splits)
    close(s)

    # Dims file
    d = open("$(prefix)_d.bin", "w+")
    write(d, size(ta.dims,1))
    write(d, size(ta.dims,2))
    write(d, ta.dims)
    close(d)

    # Leaf_data file
    l = open("$(prefix)_l.bin", "w+")
    write(l, size(ta.leaf_data,1))
    write(l, size(ta.leaf_data,2))
    write(l, ta.leaf_data)
    close(l)
end

# ḣ₁splits = [0, 30, 50]
# ḣ₁widths = [3, 5, 10]

function discretize_for_mdp(root_node::VERIFYNODE, init_lbs, init_ubs)

    # Stacks to keep track of nodes to check
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{Union{VERIFYNODE}}()

    # Add the root node to the stacks
    push!(lb_s, init_lbs)
    push!(ub_s, init_ubs)
    push!(s, root_node)

    while length(s) > 0
        #println(length(s))
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)
        curr_node = pop!(s)

        if typeof(curr_node.left) == LEAFNODE
            discretize!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs)
        else
            # Get new bounds
            # Go left, upper bounds will change
            left_lbs = copy(curr_lbs)
            left_ubs = copy(curr_ubs)
            left_ubs[dim] = split
            # Go right, lower bounds will change
            right_lbs = copy(curr_lbs)
            right_lbs[dim] = split
            right_ubs = copy(curr_ubs)

            # Add everything to the stacks
            push!(s, curr_node.left)
            push!(lb_s, curr_lbs)
            push!(ub_s, left_ubs)
            push!(s, curr_node.right) 
            push!(lb_s, right_lbs)
            push!(ub_s, curr_ubs)
        end
    end

    return root_node
end

function discretize!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs)
    ubs = unnormalize_bounds(curr_ubs)
    lbs = unnormalize_bounds(curr_lbs)
    
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
        split_specific_dims!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, dims_to_split)
    end
end

function discretize_for_mdp_3d(root_node::VERIFYNODE, init_lbs, init_ubs)

    # Stacks to keep track of nodes to check
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{Union{VERIFYNODE}}()

    # Add the root node to the stacks
    push!(lb_s, init_lbs)
    push!(ub_s, init_ubs)
    push!(s, root_node)

    while length(s) > 0
        #println(length(s))
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)
        curr_node = pop!(s)

        if typeof(curr_node.left) == LEAFNODE
            discretize_3d!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs)
        else
            # Get new bounds
            # Go left, upper bounds will change
            left_lbs = copy(curr_lbs)
            left_ubs = copy(curr_ubs)
            left_ubs[dim] = split
            # Go right, lower bounds will change
            right_lbs = copy(curr_lbs)
            right_lbs[dim] = split
            right_ubs = copy(curr_ubs)

            # Add everything to the stacks
            push!(s, curr_node.left)
            push!(lb_s, curr_lbs)
            push!(ub_s, left_ubs)
            push!(s, curr_node.right) 
            push!(lb_s, right_lbs)
            push!(ub_s, curr_ubs)
        end
    end

    return root_node
end

function discretize_3d!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs)
    ubs = unnormalize_bounds(curr_ubs)
    lbs = unnormalize_bounds(curr_lbs)
    
    # Decide if need to split
    h_inrange_pos = lbs[1] .< hsplits .< ubs[1]
    h_inrange_neg = lbs[1] .< -hsplits .< ubs[1]
    h_inrange = h_inrange_pos .| h_inrange_neg

    ḣ₀_inrange_pos = lbs[2] .< ḣ₀splits .< ubs[2]
    ḣ₀_inrange_neg = lbs[2] .< -ḣ₀splits .< ubs[2]
    ḣ₀_inrange = ḣ₀_inrange_pos .| ḣ₀_inrange_neg

    ḣ₁_inrange_pos = lbs[2] .< ḣ₁splits .< ubs[2]
    ḣ₁_inrange_neg = lbs[2] .< -ḣ₁splits .< ubs[2]
    ḣ₁_inrange = ḣ₁_inrange_pos .| ḣ₁_inrange_neg

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
    ḣ₁bin = findfirst(ḣ₀_inrange)
    if ḣ₁bin == nothing
        ḣ₁bin = findfirst(max(abs(lbs[2]), abs(ubs[2])) .< ḣ₁splits)
    end

    dims_to_split = []
    curr_diff = ubs .- lbs
    curr_diff[1] > hwidths[hbin] ? push!(dims_to_split, 1) : nothing
    curr_diff[2] > ḣ₀widths[ḣ₀bin] ? push!(dims_to_split, 2) : nothing
    curr_diff[3] > ḣ₁widths[ḣ₀bin] ? push!(dims_to_split, 3) : nothing

    next_nodes = []

    if length(dims_to_split) > 0
        split_specific_dims!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, dims_to_split)
    end
end

# """
# KD Trees
# """

# # mutable struct KDNODE
# #     dim::Int32
# #     split::Float64
# #     left::KDNODE
# #     right::KDNODE
# # end

# function equivalent_KD(curr_node::NODE)
#     children = curr_node.children
#     lbs = curr_node.lbs
#     ubs = curr_node.ubs

#     lb_diffs = [child.lbs .- lbs for child in children]

#     # Find out where the splits are
#     ndims = length(lbs)
#     dims = []
#     splits = []
#     for i = 1:length(lb_diffs)
#         diff_dims = findall(lb_diffs[i] .> 0.0)
#         for dim in diff_dims
#             if !(dim in dims)
#                 push!(dims, dim)
#                 push!(splits, children[i].lbs[dim])
#             end
#         end
#     end
# end

# # Hard-coded for 2D with splitting strategy in all dimensions
# function equivalent_2D(curr_node::NODE)
#     children = curr_node.children
#     if length(children) == 0
#         if length(curr_node.cats) == 0
#             curr_node.cats = collect(0:8)
#             get_categories!(curr_node, network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra03_v5_25HU_1000.nnet")
#         end
#         return leafnode(cats = curr_node.cats)
#     elseif length(children) == 2
#         # Figure out which dimension got split
#         lbdiffs = children[2].lbs .- curr_node.lbs
#         split_dim = findfirst(lbdiffs .!= 0.0)
#         split = children[2].lbs[split_dim]

#         # Create the tree
#         top = kdnode(dim = split_dim, split = split)

#         # Recursively add children
#         top.left = equivalent_2D(children[1])
#         top.right = equivalent_2D(children[2])
#         return top
#     else
#         first_split = children[2].lbs[1]
#         second_split = children[3].lbs[2]

#         # Create the tree
#         top = kdnode(dim = 1, split = first_split)
#         top.left = kdnode(dim = 2, split = second_split)
#         top.right = kdnode(dim = 2, split = second_split)

#         # Recursively add children
#         top.left.left = equivalent_2D(children[1])
#         top.right.left = equivalent_2D(children[2])
#         top.left.right = equivalent_2D(children[3])
#         top.right.right = equivalent_2D(children[4])
#         return top
#     end
# end

# function equivalent_2D(ver_node::VERIFYNODE)
#     #
# end

# function equivalent_2D_helper(ver_node::VERIFYNODE, curr_node::KDNODE)
#     # Base case: reached leaf node
#     if typeof(ver_node.left) == leafnode()
        
#     end

#     # Recursive case, traverse tree
#     curr_node.left = kdnode()

