using DataStructures

include("tree_utils.jl")
include("run_reluval.jl")
include("nnet_functions.jl")

function adaptive_verification(init_lbs, init_ubs; 
    min_diff = 0.1, max_dim = 2, network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    
    # Stacks to keep track of nodes to check
    # NOTE: stacks make it depth first, which I think
    # will keep the stack shorter, but can change to
    # queue if want breadth first
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{Union{VERIFYNODE}}()

    # Add the root node to the stacks
    push!(lb_s, init_lbs)
    push!(ub_s, init_ubs)
    root_node = verifynode()
    push!(s, root_node)

    while length(s) > 0
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)
        curr_node = pop!(s)

        # For now, just split in all dimensions
        verify_and_split!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs,
            min_diff = min_diff, network_path = network_path)
    end

    return root_node
end

function verify_and_split!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs; 
    min_diff = 0.1, network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    
    

    # First run verification on possible advisories
    new_cats = get_categories(curr_node.cats, curr_lbs, curr_ubs, network_path = network_path)
    curr_node.cats = new_cats
    # Decide if we should split
    max_diff = maximum(curr_ubs .- curr_lbs)
    split = (length(curr_node.cats) > 1) && (max_diff > min_diff)
    if split
        split_both!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs)
    end
end

function split!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, dim)
    split = (curr_lbs[dim] + curr_ubs[dim]) / 2

    curr_node.split = split
    curr_node.dim = dim

    curr_node.left = verifynode(cats = curr_node.cats)
    curr_node.right = verifynode(cats = curr_node.cats)

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

function split_both!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs)
    # Split in first dimension and do not add to stack
    dim = 1
    split = (curr_lbs[dim] + curr_ubs[dim]) / 2

    curr_node.split = split
    curr_node.dim = dim

    curr_node.left = verifynode(cats = curr_node.cats)
    curr_node.right = verifynode(cats = curr_node.cats)

    # Get new bounds
    # Go left, upper bounds will change
    left_lbs = copy(curr_lbs)
    left_ubs = copy(curr_ubs)
    left_ubs[dim] = split
    # Go right, lower bounds will change
    right_lbs = copy(curr_lbs)
    right_lbs[dim] = split
    right_ubs = copy(curr_ubs)

    split!(s, lb_s, ub_s, curr_node.left, left_lbs, left_ubs, 2)
    split!(s, lb_s, ub_s, curr_node.right, right_lbs, right_ubs, 2)
end

function evaluate_corners(nnet, lbs, ubs; max_dim = 2, pra = 1)
    ncorners = 2^max_dim
    points = zeros(4, ncorners)

    lbs, ubs = unnormalize_bounds(lbs), unnormalize_bounds(ubs)
    diffs = ubs .- lbs

    for i = 0:ncorners-1
        corner = digits(i, base = 2, pad = 4)
        points[:, i+1] = lbs .+ corner .* diffs
    end

    vals = evaluate_network_multiple(nnet, points)
    return [argmax(vals[:,i]) for i = 1:ncorners] .- 1
end