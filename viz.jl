#  Functions for visualization
using PGFPlots
using Colors
using ColorBrewer

include("nnet_functions.jl")
include("tree_utils.jl")

COC = 0
DNC = 1
DND = 2
DES1500 = 3
CL1500 = 4
SDES1500 = 5
SCL1500 = 6
SDES2500 = 7
SCL2500 = 8

actions = [COC, DNC, DND, DES1500, CL1500, SDES1500, SCL1500, SDES2500, SCL2500]

function cells_on_nn(root_node, depth; 
    network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")

    # Colors
    ra_1 = RGB(1.,1.,1.) # white
    ra_2 = RGB(.0,1.0,1.0) # cyan
    ra_3 = RGB(144.0/255.0,238.0/255.0,144.0/255.0) # lightgreen
    ra_4 = RGB(30.0/255.0,144.0/255.0,1.0) # dodgerblue
    ra_5 = RGB(0.0,1.0,.0) # lime
    ra_6 = RGB(0.0,0.0,1.0) # blue
    ra_7 = RGB(34.0/255.0,139.0/255.0,34.0/255.0) # forestgreen
    ra_8 = RGB(0.0,0.0,128.0/255.0) # navy
    ra_9 = RGB(0.0,100.0/255.0,0.0) # darkgreen

    colors = [ra_1;ra_2;ra_3;ra_4;ra_5;ra_6;ra_7;ra_8;ra_9]

    nnet = read_network(network_path)

    root_lbs = unnormalize_bounds(root_node.lbs)
    ḣ₁ = root_lbs[3]
    τ = root_lbs[4]

    xmin = -8000
    xmax = 8000
    ymin = -100
    ymax = 100
    nbins = 100

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="8cm", 
    xlabel=L"$h$", ylabel=L"$\dot{h}_0$", title="Neural Network Advisories")

    # Policy
    inputsNet = hcat([[h, ḣ₀, ḣ₁, τ] for h=LinRange(xmin,xmax,nbins) for ḣ₀=LinRange(ymin,ymax,nbins)]...)
    q_nnet = evaluate_network_multiple(nnet, inputsNet)

    ind = 1
    function get_heat_nn(x, y)
        qvals = q_nnet[:,ind]
        ind += 1
        return actions[argmax(qvals)]
    end

    push!(ax, Plots.Image(get_heat_nn, (xmin, xmax), (ymin, ymax), zmin = 0, zmax = 8,
    xbins = nbins, ybins = nbins, colormap = ColorMaps.RGBArrayMap(colors), colorbar=false))

    leaves = [];
    get_leaves(root_node, leaves, depth)

    # Rectangles
    for leaf in leaves
        lbs = unnormalize_bounds(leaf.lbs)
        ubs = unnormalize_bounds(leaf.ubs)
        push!(ax, Plots.Command(get_rectangle(lbs, ubs)))
    end

    return ax
end

function get_rectangle(lb, ub)
    return "\\draw ($(string(lb[1])),$(string(lb[2]))) rectangle ($(string(ub[1])),$(string(ub[2])));"
end

function unnormalize_bounds(bounds; ranges = [16000.0, 200.0, 200.0, 40.0], means = [0.0, 0.0, 0.0, 20.0])
    return bounds .* ranges + means
end