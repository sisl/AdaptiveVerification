# for pra in pras
#     τ = 0
#     println("pra: $pra, τ: 0")

#     trees_rc = Dict()


#     root_node = verifynode()
#     init_lbs = [-0.5, -0.5, -0.5, -0.5]
#     init_ubs = [0.5, 0.5, 0.5, -0.5]
#     discretize_for_mdp_3d(root_node, init_lbs, init_ubs)
#     ta = write_to_arrays(root_node)
#     tas[(pra - 1, 0.0)] = ta

#     # for normalized_τ in normalized_τs
#     #     τ = normalized_τ * 40 + 20
#     #     println("pra: $pra, τ: $τ")
#     #     start = time()
#     #     tree = kdtrees[(pra - 1, τ)]
#     #     ta = write_to_arrays(tree)
#     #     tas[(pra - 1, τ)] = ta
#     #     println("Time required: $(time() - start)")
#     # end
# end

# trees_rc = Dict()
# trees_final = Dict()

# for pra in 1:9
#     root_node = verifynode()
#     init_lbs = [-0.5, -0.5, -0.5, -0.5]
#     init_ubs = [0.5, 0.5, 0.5, -0.5]
#     trees_rc[pra] = remotecall(discretize_for_mdp_3d, mod(pra + 1, nprocs() - 1) + 2, root_node, init_lbs, init_ubs)
# end

# for pra in 1:9
#     trees_final[pra] = fetch(trees_rc[pra])
# end

# for pra in 1:9
#     ta = write_to_arrays(trees_final[pra], max_dim = 3)
#     write_to_files(ta, "/scratch/smkatz/pre3d/pra$(pra)tau0")
# end

τs = collect(1:1)
normalized_τs = (τs .- 20) ./ 40
pras = collect(1:9)
min_diff = 0.3

for normalized_τ in normalized_τs
    τ = convert(Int64, normalized_τ * 40 + 20)
    
    trees_rc = Dict()
    trees_final = Dict()

    init_lbs = [-0.5, -0.5, -0.5, normalized_τ]
    init_ubs = [0.5, 0.5, 0.5, normalized_τ]

    for pra in pras
        network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0$(pra)_v5_25HU_1000.nnet"
        trees_rc[pra] = remotecall(adaptive_verification, mod(pra + 1, nprocs() - 1) + 2, init_lbs, init_ubs,
                            min_diff = min_diff, max_dim = 3, network_path = network_path)
    end

    for pra in pras
        trees_final = fetch(trees_rc[pra])
    end

    for pra in pras
        ta = write_to_arrays(trees_final[pra], max_dim = 3)
        write_to_files(ta, "/scratch/smkatz/pre3d/pra$(pra)tau$(τ)")
    end
end
