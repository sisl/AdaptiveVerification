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