include("opt1.jl")
include("ss2.jl")
using Ipopt

ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer)

function update_pipe_fields!(ref)
    for (i, pipe) in ref[:nw][0][:pipe]
        pd_min = pipe["pd_min"]
        pd_max = pipe["pd_max"]
        lambda = pipe["friction_factor"]
        L = pipe["length"] * ref[:base_length]
        D = pipe["diameter"]
        pipe["resistance"] = lambda * L / D
        w = 1 / pipe["resistance"]
        min_flux = pd_min < 0 ? -sqrt(w * abs(pd_min)) : sqrt(w * abs(pd_min))
        max_flux = pd_max < 0 ? -sqrt(w * abs(pd_max)) : sqrt(w * abs(pd_max))
        pipe["flux_min"] = min_flux
        pipe["flux_max"] = max_flux
    end

    for (i, pipe) in ref[:nw][0][:compressor]
        pd_min = pipe["pd_min"]
        pd_max = pipe["pd_max"]
        lambda = pipe["friction_factor"]
        L = pipe["length"] * ref[:base_length]
        D = pipe["diameter"]
        pipe["resistance"] = lambda * L / D
        w = 1 / pipe["resistance"]
        min_flux = pd_min < 0 ? -sqrt(w * abs(pd_min)) : sqrt(w * abs(pd_min))
        max_flux = pd_max < 0 ? -sqrt(w * abs(pd_max)) : sqrt(w * abs(pd_max))
        pipe["flux_min"] = min_flux
        pipe["flux_max"] = max_flux
        if pipe["directionality"] == 1
            pipe["flux_min"] = 0.0
            min_flux = 0.0
        end
        pipe["flow_min"] = pipe["area"] * min_flux
        pipe["flow_max"] = pipe["area"] * max_flux
    end

    slack_junctions = Dict()
    for (k, v) in ref[:nw][0][:junction]
        (v["junction_type"] == 1) && (slack_junctions[k] = v)
    end

    ref[:nw][0][:slack_junctions] = slack_junctions

end

struct SSModel
    data
    ref
    model
    var
    con
    sol
end

function main2()
    data = GasModels.parse_file("../data/case6a.m")
    ref = GasModels.build_ref(data)
    update_pipe_fields!(ref)
    partition =[[1,5],[2,3,4,6]] # these divide the nodes in the network
    # based upon the partition they belong to. currently the number of partitions
    # is 2 each having 2 and 4 nodes (junctions) respectively. Note p1∪ p2 = 𝐕
    # the common components are pipes, compressors never junctions. the dec variables
    # which are common to both networks will be pressure, flow alpha and phi. the supply
    # and demand points. m, var, con, expr, component = post_ss_model(ref,partition)
    β = 1000.0
    λ_h= 1000000
    λ_l = -1000000
    γ= 2.0
    ρ = γ*β
    ω= 0.75
    k =1
    t =1
    C=C1=C2=C3=100
    ϵ=0.00001
    ϵ3= -1.0
    model = post_ss_model(ref, partition, β)
    modZ0= model[4][1][:params][:outertolerance]
    #proj = [λ_h- λ_l for i in 1: length(λ)]
    l = length(partition)
    #ϵ1= ϵ2= ϵ3= 2*sqrt(q)/(k*ρ)
    print("SET MODEL DONE \n")
    while(modZ0>ϵ)
        while(C3> ϵ3)
            print("ENTER INNER LOOP ", t, "\n")
            if (t>1 ||k>1)
                optOver!("X", model)
            end
            s= "local"
            model = L(model, ref, ρ,β, s) # optimizing over X-local
            print("Ater OPT 1 ",k,"\t", t, "\n")
            optOver!("Global", model)
            s="global"
            model = L(model, ref, ρ,β, s)# optimizing over X-global
            print("Ater OPT 2 ",k,"\t", t, "\n")
            optOver!("Z", model)
            s = "z"
            model = L(model, ref, ρ,β,s)# optimizing over z variables
            print("Ater OPT 3 ",k,"\t", t, "\n")
            updateY!(model,ref,ρ, k) # dual update of Y
            print("AfterY update ", t, "\n")
            C3= model[4][1][:params][:tolerance]
            print("inner tolerance \t", C3, "\n")
            ϵ3 = model[4][1][:params][:tolerancelimit]
            print("inner tolerancelimit \t", ϵ3,"\n")
            t+=1
            print(t,"\n")
            # if t>2
            #     return
            # end
        end
        updateλ!(model,β, λ_h, λ_l)
        modZ= model[4][1][:params][:outertolerance]
        print("outer tolerance\t", modZ, "\n")
         if modZ>= modZ0*ω
             β = γ*β
         else
             β = γ*β
         end
         modZ0= modZ

        C3 = 100
        ϵ3 =-1
        ρ = γ*β
        k +=1
        t=1
        reinitialize!(model, β,ref)
        print(k)
         # if(k==3)
         #     return
         # end
    end
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n")
    sol = Dict()
    vname= ["junction","pipe","receipt","delivery","compressor","compressor"]
    vname2 = ["pressure","pipe_flux","fs","fd","compressor_flow","compressor_ratio"]
    vname3 = ["p","phi","s","d","f","α"]
    sol[:pressure]= [Dict() for z in 1:l+1]
    sol[:fs]= [Dict() for z in 1:l+1]
    sol[:fd]= [Dict() for z in 1:l+1]
    for z in 1:l+1
        for (k,v) in model[4][z][:Xp]
            sol[:pressure][z][k] = v*  ref[:base_pressure]
        end
        for (k,v) in model[4][z][:Xs]
            sol[:fs][z][k] = v
        end
        for (k,v) in model[4][z][:Xd]
            sol[:fd][z][k] = v
        end
    end

    # for (i, pipe) in ref[:nw][0][:pipe]
    #     sol[:pipe_flow][i] = model[4][l+1][:Xphi][i] * pipe["area"] * ref[:base_flux]
    # end
    #
    # for (i, compressor) in ref[:nw][0][:compressor]
    #     fr = compressor["fr_junction"]
    #     to = compressor["to_junction"]
    #     sol[:compressor_flow][i] = model[4][l+1][:Xf][i] * ref[:base_flow]
    #     sol[:compressor_ratio][i] = model[4][l+1][:Xα][i]
    #     #sol[:compressor_power][i] = model[4][l+1][:compressor_power][i] * ref[:base_flow]
    #     sol[:compressor_inlet_pressure][i] = sol[:pressure][fr]
    #     sol[:compressor_discharge_pressure][i] = sol[:pressure][to]
    # end
    #
    # for (i, receipt) in ref[:nw][0][:receipt]
    #     sol[:supply][i] = model[4][l+1][:Xs][i] * ref[:base_flow]
    # end
    #
    # for (i, delivery) in ref[:nw][0][:delivery]
    #     sol[:demand][i] = model[4][l+1][:Xd][i] * ref[:base_flow]
    # end
    for z in 1:l+1
        for (k,v) in sol[:pressure][z]
            if z<=l
                print("local ", z, "\t", k,"\t",v, "\n")
            else
                print("global ", z, "\t", k,"\t",v, "\n")
            end

        end
        for (k,v) in sol[:fs][z]
            if z<=l
                print("local ", z, "\t", k,"\t",v, "\n")
            else
                print("global ", z, "\t", k,"\t",v, "\n")
            end

        end
        for (k,v) in sol[:fd][z]
            if z<=l
                print("local ", z, "\t", k,"\t",v, "\n")
            else
                print("global ", z, "\t", k,"\t",v, "\n")
            end

        end
    end
    print("opt ob val =  ", model[4][1][:params][:objVal], "\n")
    return SSModel(data, ref, model[1], var, model[3], sol)
end
