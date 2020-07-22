function get_categories!(curr_node::NODE; network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    lbs, ubs = unnormalize_bounds(curr_node.lbs), unnormalize_bounds(curr_node.ubs)
    old_cats = deepcopy(curr_node.cats)
    times, poss_cats = verify_region(lbs, ubs, advs = old_cats, network_path = network_path)
    poss_inds = findall(poss_cats .> 0)
    curr_node.cats = old_cats[poss_inds]
end

function unnormalize_bounds(bounds; ranges = [16000.0, 200.0, 200.0, 40.0], means = [0.0, 0.0, 0.0, 20.0])
    return bounds .* ranges + means
end

function verify_region(lbs, ubs; advs = collect(0:8), network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    poss_advs = []
    advs_string = "$(advs[1])"
    nadvs = length(advs)
    for i = 2:nadvs
        advs_string = "$advs_string $(advs[i])"
    end
    cmd = `./network_test_VC $(network_path) $(lbs[1]) $(lbs[2]) $(lbs[3]) $(lbs[4]) $(ubs[1]) $(ubs[2]) $(ubs[3]) $(ubs[4]) $advs`
    output = read(cmd, String)
    vals = parse.(Float64, split(output, ",")[1:end-1])
    times = vals[2:2:2*nadvs]
    advs = vals[1:2:2*nadvs-1]
    return times, advs
end