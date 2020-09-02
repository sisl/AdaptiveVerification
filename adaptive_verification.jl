using DataStructures

include("tree_utils.jl")
include("run_reluval.jl")
include("nnet_functions.jl")

function adaptive_verification(init_lbs, init_ubs; 
    min_diff = 0.1, max_dim = 2, network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    
    nnet = read_network(network_path)

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
    # Run verification on the initial node
    new_cats = get_categories(root_node.cats, init_lbs, init_ubs, network_path = network_path)
    root_node.cats = new_cats
    push!(s, root_node)

    index = 1

    while length(s) > 0
        #println(length(s))
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)
        curr_node = pop!(s)

        verify_and_split!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, nnet,
            max_dim = max_dim, min_diff = min_diff, network_path = network_path)
    end

    return root_node
end

function verify_and_split_mdp!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, nnet; 
    min_diff = 0.1, max_dim = 2, network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    
    # Evaluate the corners
    corner_advs = evaluate_corners_faster(nnet, curr_lbs, curr_ubs, max_dim = max_dim)
    # Decide if we should split
    max_diff = minimum(curr_ubs[1:max_dim] .- curr_lbs[1:max_dim])
    if (length(unique(corner_advs)) > 1) && (max_diff > min_diff)
        split_dims = split_dims_from_corners_3d(corner_advs)
        split_specific_dims!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, split_dims)
    else
        # Run verification on possible advisories
        new_cats = get_categories(curr_node.cats, curr_lbs, curr_ubs, network_path = network_path)
        curr_node.cats = new_cats
        split = (length(curr_node.cats) > 1) && (max_diff > min_diff)
        if split
            split_multiple!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, max_dim)
        end
    end
end

function verify_and_split!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, nnet; 
    min_diff = 0.1, max_dim = 2, network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    
    # Evaluate the corners
    corner_advs = evaluate_corners_faster(nnet, curr_lbs, curr_ubs, max_dim = max_dim)
    # Decide if we should split
    max_diff = minimum(curr_ubs[1:max_dim] .- curr_lbs[1:max_dim])
    if (length(unique(corner_advs)) > 1) && (max_diff > min_diff)
        split_dims = split_dims_from_corners_3d(corner_advs)
        split_specific_dims!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, split_dims)
    else
        # Run verification on possible advisories
        new_cats = get_categories(curr_node.cats, curr_lbs, curr_ubs, network_path = network_path)
        curr_node.cats = new_cats
        split = (length(curr_node.cats) > 1) && (max_diff > min_diff)
        if split
            split_multiple!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, max_dim)
        end
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

function split_multiple!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, max_dim)
    split!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, 1)
    for i = 2:max_dim
        nodes = []
        lbs = []
        ubs = []
        for j = 1:2^(i-1)
            push!(nodes, pop!(s))
            push!(lbs, pop!(lb_s))
            push!(ubs, pop!(ub_s))
        end
        for j = 1:2^(i-1)
            split!(s, lb_s, ub_s, nodes[j], lbs[j], ubs[j], i)
        end
    end
end

function split_specific_dims!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, dims)
    split!(s, lb_s, ub_s, curr_node, curr_lbs, curr_ubs, dims[1])
    for i = 2:length(dims)
        nodes = []
        lbs = []
        ubs = []
        for j = 1:2^(i-1)
            push!(nodes, pop!(s))
            push!(lbs, pop!(lb_s))
            push!(ubs, pop!(ub_s))
        end
        for j = 1:2^(i-1)
            split!(s, lb_s, ub_s, nodes[j], lbs[j], ubs[j], dims[i])
        end
    end
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

function evaluate_corners_faster(nnet, lbs, ubs; max_dim = 2, pra = 1)
    ncorners = 2^max_dim
    points = zeros(4, ncorners)

    diffs = ubs .- lbs

    for i = 0:ncorners-1
        corner = digits(i, base = 2, pad = 4)
        points[:, i+1] = lbs .+ corner .* diffs
    end

    vals = evaluate_network_faster(nnet, points)
    return [argmax(vals[:,i]) for i = 1:ncorners] .- 1
end

function split_dims_from_corners_2d(advs)
    split_dims = []

    if any(abs.(advs[[1, 3]] - advs[[2, 4]]) .> 0)
        push!(split_dims, 1)
    end

    if any(abs.(advs[[1, 2]] - advs[[3, 4]]) .> 0)
        push!(split_dims, 2)
    end

    return split_dims
end

function split_dims_from_corners_3d(advs)
    split_dims = []

    if any(abs.(advs[[1, 3, 5, 7]] - advs[[2, 4, 6, 8]]) .> 0)
        push!(split_dims, 1)
    end

    if any(abs.(advs[[1, 2, 5, 6]] - advs[[3, 4, 7, 8]]) .> 0)
        push!(split_dims, 2)
    end

    if any(abs.(advs[[1, 2, 3, 4]] - advs[[5, 6, 7, 8]]) .> 0)
        push!(split_dims, 3)
    end

    return split_dims
end
