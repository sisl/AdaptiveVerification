"""
KDTrees
"""

mutable struct LEAFNODE
    cats::Vector{Int64}
    qvals::Vector{Float64}
end

mutable struct VERIFYNODE
    dim::Int32
    split::Float64
    cats::Vector{Int64}
    left::Union{LEAFNODE, VERIFYNODE}
    right::Union{LEAFNODE, VERIFYNODE}
end

mutable struct KDNODE
    dim::Int32
    split::Float64
    left::Union{KDNODE, LEAFNODE, VERIFYNODE}
    right::Union{KDNODE, LEAFNODE, VERIFYNODE}
end

function verifynode(;dim = 1, split = 0.0, cats = collect(0:8), left = leafnode(), right = leafnode())
    return VERIFYNODE(dim, split, cats, left, right)
end

function leafnode(;cats = collect(0:8), qvals = Vector{Float64}())
    return LEAFNODE(cats, qvals)
end

function kdnode(;dim = 1, split = 0.0, left = leafnode(), right = leafnode())
    return KDNODE(dim, split, left, right)
end

function get_leaf(root_node, point)
    curr_node = root_node
    while typeof(curr_node) != LEAFNODE
        val = point[curr_node.dim]
        split = curr_node.split
        curr_node = val < split ? curr_node.left : curr_node.right
    end
    return curr_node
end

function get_leaves(root_node)
    leaves = []
    get_leaves_helper(root_node, leaves)
    return leaves
end

function get_leaves_helper(root_node, leaves)
    if typeof(root_node.left) == LEAFNODE
        push!(leaves, root_node)
    else
        get_leaves_helper(root_node.left, leaves)
        get_leaves_helper(root_node.right, leaves)
    end
end

function get_bounds_and_cats(root_node)
    cats = []
    lbs = []
    ubs = []
    
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{VERIFYNODE}()

    push!(lb_s, [-0.5, -0.5])
    push!(ub_s, [0.5, 0.5])
    push!(s, root_node)

    while !isempty(s)
        curr = pop!(s)
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)

        if typeof(curr.left) == LEAFNODE
            push!(cats, curr.cats)
            push!(lbs, curr_lbs)
            push!(ubs, curr_ubs)
        else
            # Traverse tree and keep track of bounds
            dim = curr.dim
            split = curr.split
            # Go left, upper bounds will change
            left_ubs = copy(curr_ubs)
            left_ubs[dim] = split

            push!(lb_s, curr_lbs)
            push!(ub_s, left_ubs)
            push!(s, curr.left)

            # Go right, lower bounds will change
            right_lbs = copy(curr_lbs)
            right_lbs[dim] = split
            
            push!(lb_s, right_lbs)
            push!(ub_s, curr_ubs)
            push!(s, curr.right)
        end
    end

    return lbs, ubs, cats
end