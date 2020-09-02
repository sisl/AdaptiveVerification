#  Functions for visualization
using DataStructures
using PGFPlots
using Colors
using ColorBrewer
using Reel
Reel.extension(m::MIME"image/svg+xml") = "svg"
Reel.set_output_type("gif") # may be necessary for use in IJulia

#include("nnet_functions.jl")
#include("tree_utils.jl")

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

function cells_on_nn(lower_bound_list::Vector, upper_bound_list::Vector; 
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

    # Unnormalize everying
    lbs = [lower_bound_list[i] .* nnet.ranges[1:4] .+ nnet.means[1:4] for i = 1:length(lower_bound_list)]
    ubs = [upper_bound_list[i] .* nnet.ranges[1:4] .+ nnet.means[1:4] for i = 1:length(upper_bound_list)]

    ḣ₁ = lbs[1][3]
    τ = lbs[1][4]

    xmin = -8000
    xmax = 8000
    ymin = -100
    ymax = 100
    nbins = 100

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="7cm", height="8cm", 
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

    #Rectangles
    for i = 1:length(lbs)
        push!(ax, Plots.Command(get_rectangle(lbs[i], ubs[i])))
    end

    return ax
end

function plot_nadvs(tas::TREE_ARRAYS)
    
    lbs_u, ubs_u, cats = get_bounds_and_cats(tas)

    # Unnormalize everying
    lbs = [unnormalize_point(lbs_u[i]) for i = 1:length(lbs_u)]
    ubs = [unnormalize_point(ubs_u[i]) for i = 1:length(ubs_u)]

    ymin = -8000
    ymax = 8000
    xmin = -100
    xmax = 100

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$\dot{h}_0$", ylabel=L"$h$", title="Number of Possible Advisories")

    for i = 1:length(lbs)
        if length(cats[i]) == 1
            color = "blue"
        elseif length(cats[i]) == 2
            color = "red"
        else
            color = "yellow"
        end
        push!(ax, Plots.Command(get_filled_rectangle([lbs[i][2], lbs[i][1]],
                                                     [ubs[i][2], ubs[i][1]], color)))
    end

    return ax
end

function plot_nadvs(tree::VERIFYNODE)
    
    lbs_u, ubs_u, cats = get_bounds_and_cats(tree)

    # Unnormalize everying
    lbs = [unnormalize_point(lbs_u[i]) for i = 1:length(lbs_u)]
    ubs = [unnormalize_point(ubs_u[i]) for i = 1:length(ubs_u)]

    ymin = -8000
    ymax = 8000
    xmin = -100
    xmax = 100

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$\dot{h}_0$", ylabel=L"$h$", title="Number of Possible Advisories")

    for i = 1:length(lbs)
        if length(cats[i]) == 1
            color = "blue"
        elseif length(cats[i]) == 2
            color = "red"
        else
            color = "yellow"
        end
        push!(ax, Plots.Command(get_filled_rectangle([lbs[i][2], lbs[i][1]],
                                                     [ubs[i][2], ubs[i][1]], color)))
    end

    return ax
end

function plot_nadvs(root_node::KDNODE)

    # Colors
    one_adv = RGB(0.0,0.0,1.0) # blue
    two_adv = RGB(1.0,0.0,0.0) # red
    three_adv = RGB(1.0,1.0,0.0) # yellow

    colors = [one_adv, two_adv, three_adv]

    xmin = -8000
    xmax = 8000
    ymin = -100
    ymax = 100
    nbins = 100

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="7cm", height="8cm", 
    xlabel=L"$h$", ylabel=L"$\dot{h}_0$", title="Neural Network Advisories")

    function get_heat(x, y)
        h = x / 16000
        ḣ₀ = y / 200
        leaf = get_leaf(root_node, [h, ḣ₀])
        return length(leaf.cats)
    end

    push!(ax, Plots.Image(get_heat, (xmin, xmax), (ymin, ymax), zmin = 1, zmax = 3,
    xbins = nbins, ybins = nbins, colormap = ColorMaps.RGBArrayMap(colors), colorbar=false))

    return ax
end

function get_key()
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
    bg_colors = [RGB(1.0,1.0,1.0)]

    # Create scatter plot classes for color key
    sc_string = "{"
    for i=1:9
        define_color("ra_$i",  colors[i])
        if i==1
            sc_string *= "ra_$i={mark=square, style={black, mark options={fill=ra_$i}, mark size=6}},"
        else
            sc_string *= "ra_$i={style={ra_$i, mark size=6}},"
        end
    end

    # Color key as a scatter plot
    sc_string=sc_string[1:end-1]*"}"
    xx = [-1.5,-1.5,-1.5, -1.5, -1.5, -1.5, -1.5, 0.4 ,.4,]
    yy = [1.65,1.15,0.65, 0.15, -0.35, -0.85, -1.35, 1.65, 1.15]
    zz = ["ra_1","ra_2","ra_3","ra_4","ra_5","ra_6","ra_7","ra_8","ra_9"]
    sc = string(sc_string)

    f = (x,y)->x # Dummy function for background white image
    key = Axis([
        Plots.Image(f, (-2,2), (-2,2),colormap = ColorMaps.RGBArrayMap(bg_colors),colorbar=false),
        Plots.Scatter(xx, yy, zz, scatterClasses=sc),
        Plots.Node("RA 1: COC ",0.15,0.915,style="black,anchor=west", axis="axis description cs"),
        Plots.Node("RA 2: DNC ",0.15,0.790,style="black,anchor=west", axis="axis description cs"),
        Plots.Node("RA 3: DND",0.15,0.665,style="black,anchor=west", axis="axis description cs"),
        Plots.Node("RA 4: DES15000",0.15,0.540,style="black,anchor=west", axis="axis description cs"),
        Plots.Node("RA 5: CL1500 ",0.15,0.415,style="black,anchor=west", axis="axis description cs"),
        Plots.Node("RA 6: SDES1500",0.15,0.290,style="black,anchor=west", axis="axis description cs"),
        Plots.Node("RA 7: SCL1500",0.15,0.165,style="black,anchor=west", axis="axis description cs"),
        Plots.Node("RA 8:  SDES2500",0.63,0.915,style="black,anchor=west", axis="axis description cs"),
        Plots.Node("RA 9:  SCL2500",0.63,0.790,style="black,anchor=west", axis="axis description cs"),
        ],width="10cm",height="8cm", hideAxis =true, title="KEY")

    return key
end

function plot_all(lower_bound_list::Vector, upper_bound_list::Vector, cats;
    network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    
    ax1 = plot_nadvs(lower_bound_list, upper_bound_list, cats, network_path = network_path)
    ax2 = cells_on_nn(lower_bound_list, upper_bound_list, network_path = network_path)
    key = get_key()
    
    g = GroupPlot(3, 1, groupStyle = "horizontal sep=2cm")
    push!(g, ax1)
    push!(g, ax2)
    push!(g, key)
    g
end

""" Creating gifs """

""" Plotting utils """

function get_rectangle(lb, ub)
    return "\\draw ($(string(lb[1])),$(string(lb[2]))) rectangle ($(string(ub[1])),$(string(ub[2])));"
end

function get_filled_rectangle(lb, ub, color)
    return "\\filldraw[fill=$(color), draw=black] ($(string(lb[1])),$(string(lb[2]))) rectangle ($(string(ub[1])),$(string(ub[2])));"
end

function unnormalize_bounds(bounds; ranges = [16000.0, 200.0, 200.0, 40.0], means = [0.0, 0.0, 0.0, 20.0])
    return bounds .* ranges + means
end

function unnormalize_point(point::Vector{Float64})
    return point .* [16000.0, 200.0]
end