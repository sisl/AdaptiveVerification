
mutable struct NODE
    lbs::Vector{Float64}
    ubs::Vector{Float64}
    cats::Vector{Int64}
    children::Vector{NODE}
end

function node(;lbs = Vector{Float64}(), 
    ubs = Vector{Float64}(), 
    cats = collect(0:8), 
    children = Vector{NODE}())
    return NODE(lbs, ubs, cats, children)
end

# Modifies leaves (not the tree)
function get_leaves(curr_node, leaves, depth)
    if depth == 0
        push!(leaves, curr_node)
        return
    else
        if length(curr_node.children) > 0
            for child in curr_node.children
                get_leaves(child, leaves, depth - 1)
            end
        else
            push!(leaves, curr_node)
        end
    end
end

function max_depth(curr_node)
    if length(curr_node.children) == 0
        return 1
    else
        depths = zeros(length(curr_node.children))
        for i = 1:length(depths)
            depths[i] = max_depth(curr_node.children[i])
        end
        return 1 + maximum(depths)
    end
end

function get_bounds_list(root_node)
    q = Queue{NODE}()
    enqueue!(q, root_node)
    lbs = []
    ubs = []
    cats = []
    while length(q) > 0
        curr_node = dequeue!(q)
        push!(lbs, curr_node.lbs)
        push!(ubs, curr_node.ubs)
        push!(cats, curr_node.cats)
        children = curr_node.children
        for child in children
            enqueue!(q, child)
        end
    end
    return lbs, ubs, cats
end