include("part.jl")
include("set.jl")
using JuMP
using GasModels
# this module addds all the variables and constraints for the local problems and
# their objectives to the model once and for all so we can apply the algorithm to it
function post_ss_model(ref, partition, β)
    nw = ref[:nw][0]
    l = length(partition)
    part = partition1(partition, ref)
    m,var, con, expr, arr = setVariables(nw, l, part, β)
    p,supply,demand,phi,flow,alpha,Yp,Ysupply,Ydemand,Yphi,Yflow,
    Yalpha, Zp,Zsupply,Zdemand,Zphi,Zflow,Zalpha= arr#,Sp, Sphi, Sflow, Salpha= arr
    # creating power expression local variable

    for z in 1:l
         ct=0
         for i in part[1][z][:compressor]
                a = alpha[z][i]#i in part[2][z][:compressor] ? Salpha[z][i] : alpha[z][i]
                f = flow[z][i]#i in part[2][z][:compressor] ? Sflow[z][i] : flow[z][i]
                m1 = (ref[:specific_heat_capacity_ratio] - 1) / ref[:specific_heat_capacity_ratio]
                W = 286.76 * ref[:temperature] / ref[:gas_specific_gravity] / m1
                var[z][:compressor_power][i] = JuMP.@NLexpression(m, W * abs(f) * (a^m1 - 1.0))
                ct+=1
        end
        print("Local power variables created in partition ",z, " is\t ", ct, "\n")
    end

    print("power expressions created \n")

    for z in 1:l
        ct=0
        for i in partition[z]
            var[z][:net_nodal_edge_out_flow][i] = 0
                ct+=1
                for j in nw[:pipes_fr][i]
                    #flux = j in part[2][z][:pipe] ? Sphi[z][j] : phi[z][j]
                    var[z][:net_nodal_edge_out_flow][i] += (phi[z][j] * nw[:pipe][j]["area"])
                end
                for j in nw[:compressors_fr][i]
                    #fl = j in part[2][z][:compressor] ? Sflow[z][j] : flow[z][j]
                    var[z][:net_nodal_edge_out_flow][i] += flow[z][j]
                end
                for j in nw[:pipes_to][i]
                    #flux = j in part[2][z][:pipe] ? Sphi[z][j] : phi[z][j]
                    var[z][:net_nodal_edge_out_flow][i] -= (phi[z][j] * nw[:pipe][j]["area"])
                end
                for j in nw[:compressors_to][i]
                    #fl = j in part[2][z][:compressor] ? Sflow[z][j] : flow[z][j]
                    var[z][:net_nodal_edge_out_flow][i] -= flow[z][j]
                end
        end
        print("Local flow nodal outflow variables created in partition ",z, " is\t ", ct, "\n")
    end


    for z in 1:l
        ct=0
        for i in partition[z]
            var[z][:net_nodal_injection][i] =0
            ct+=1
            for j in nw[:dispatchable_receipts_in_junction][i]
                var[z][:net_nodal_injection][i] += supply[z][j]
            end
            for j in nw[:dispatchable_deliveries_in_junction][i]
                var[z][:net_nodal_injection][i] -= demand[z][j]
            end
            for j in nw[:nondispatchable_receipts_in_junction][i]
                var[z][:net_nodal_injection][i] += nw[:receipt][j]["injection_nominal"]
            end
            for j in nw[:nondispatchable_deliveries_in_junction][i]
                var[z][:net_nodal_injection][i] -= nw[:delivery][j]["withdrawal_nominal"]
            end
        end
        print("Local flow injection variables created in partition ",z, " is\t ", ct, "\n")
    end



    print("flow expressions created\n")
    # slack pressure on a local node
    slack = keys(nw[:slack_junctions])
    slackjunctions = [nw[:slack_junctions][s]["id"] for s in slack]

    for z in 1 :l
        ct=0
        for s in slackjunctions
            if s in part[1][z][:junction]
                con[z][:slack_pressure][s] = JuMP.@constraint(m, p[z][s] == nw[:slack_junctions][s]["p_nominal"])
                ct+=1
            end
        end
        print(ct,"\t","slack pressure constraints set in partition ", z, "\n")
    end


    # nodal balance
    for z in 1:l
        ct=0
        for i in partition[z]
            net_injection = var[z][:net_nodal_injection][i]
            net_nodal_edge_out_flow = var[z][:net_nodal_edge_out_flow][i]
            con[z][:nodal_balance][i] =JuMP.@constraint(m, net_injection == net_nodal_edge_out_flow)
            ct+=1
            # net_nodal_edge_out_flow = var[z][:net_nodal_edge_out_flowE][i]
            # con[z][:nodal_balanceE][i] =JuMP.@constraint(m, net_injection == net_nodal_edge_out_flow)
            # ct+=1
        end
        print(ct, "\t local nodal balance constraints set in partition ",z, "\n")
    end
    ct=0

    print("nodal balance set\n")

    for z in 1 :l
        ct=0
        for i in part[1][z][:pipe]
            from = nw[:pipe][i]["fr_junction"]
            to = nw[:pipe][i]["to_junction"]
            p_fr = p[z][from]#from in partition[z] ? p[z][from] : Sp[z][from]
            p_to = p[z][to]#to in partition[z] ? p[z][to] : Sp[z][to]
            resistance = nw[:pipe][i]["resistance"]
            f = phi[z][i]#i in part[2][z][:pipe] ? Sphi[z][i] : phi[z][i]
            con[z][:pipe_physics][i] =JuMP.@NLconstraint(m, p_fr^2 - p_to^2 - resistance * f * abs(f) == 0)
            ct+=1
        end
        print(ct, "\t local shared pipe physics constraints set in partition ",z, "\n")
    end

    print("pipe physics set\n")
    #compressor constraints, boost flow and max
    for z in 1:l
        ct=0
        for i in part[1][z][:compressor]
            fr = nw[:compressor][i]["fr_junction"]
            to = nw[:compressor][i]["to_junction"]
            p_fr = p[z][fr]#fr in partition[z] ? p[z][fr] : Sp[z][fr]
            p_to = p[z][to]#to in partition[z] ? p[z][to] : Sp[z][to]
            a = alpha[z][i]#i in part[2][z][:compressor] ? Salpha[z][i] : alpha[z][i]
            f = flow[z][i]#i in part[2][z][:compressor] ? Sflow[z][i] : flow[z][i]
            power = var[z][:compressor_power][i]
            power_max = nw[:compressor][i]["power_max"]
            con[z][:compressor_boost][i] = JuMP.@constraint(m, p_to == a * p_fr)
            ct+=1
            con[z][:compressor_flow_dir][i] = JuMP.@constraint(m, f * (p_fr - p_to) <= 0)
            ct+=1
            con[z][:compressor_power][i] = JuMP.@NLconstraint(m, power <= power_max)
            ct+=1
        end
        print(ct, "\t local compressor constraints in partition ",z, "\n")
    end
    print("compressor constraints set\n")
    # equality constraints
    # for z in 1:l
    #                 for j in part[2][z][:junction]
    #                         con[z][:equalityP][j]=@constraint(m, var[z][:pressure][j] ==var[l+1][:pressure][j])
    #                 end
    #                 for j in  part[2][z][:pipe]
    #                          con[z][:equalityPhi][j] =  @constraint(m, var[z][:pipe_flux][j] ==var[l+1][:pipe_flux][j])
    #                 end
    #                 for j in part[2][z][:compressor]
    #                          con[z][:equalityAlpha][j] =  @constraint(m, var[z][:compressor_ratio][j] ==var[l+1][:compressor_ratio][j])
    #                          con[z][:equalityFlow][j] =  @constraint(m, var[z][:compressor_flow][j] ==var[l+1][:compressor_flow][j])
    #                 end
                    # for j in  part[3][:delivery]#nw[:delivery]
                    #         if j in part[1][z][:delivery]
                    #                 con[z][:equalityD][j] =  @constraint(m, var[z][:fd][j] ==var[l+1][:fd][j])
                    #         end
                    # end
                    #  for j in  part[3][:receipt]#nw[:receipt]
                    #          if j in part[1][z][:receipt]
                    #                  con[z][:equalityS][j] =  @constraint(m, var[z][:fs][j] ==var[l+1][:fs][j])
                    #          end
                    # end
    #end
    #
    # create f(x) and the other objective function expressions
    #f(x) modify if to double count.. for eg use part[1][z] instead of partition to double count. partition
    # and  use of compressor for junction condition coiunts the contribution to the obj once...
    econ_weight = ref[:economic_weighting]
    load_shed_expressions = []
    compressor_power_expressions = []
    for z in 1:l
        for (i, receipt) in nw[:dispatchable_receipt]
            if (receipt["junction_id"] in partition[z] )# part[1][z][:junction] &&!(receipt["junction_id"] in part[2][z][:junction]))
                push!(load_shed_expressions, JuMP.@NLexpression(m, receipt["offer_price"] * supply[z][i]))
            end
        end
        for (i, delivery) in nw[:dispatchable_delivery]
            if (delivery["junction_id"] in partition[z])#in part[1][z][:junction] &&!(delivery["junction_id"] in part[2][z][:junction]))
                push!(load_shed_expressions,JuMP.@NLexpression(m, -delivery["bid_price"] * demand[z][i]))
            end
        end
        for i in part[1][z][:compressor]
            if (nw[:compressor][i]["fr_junction"] in partition[z])
                #&&!(nw[:compressor][i]["fr_junction"] in part[2][z][:junction]))
                push!(compressor_power_expressions, JuMP.@NLexpression(m,var[z][:compressor_power][i]))
            end
        end
    end
    expr[1][:params][:load]= load_shed_expressions
    expr[1][:params][:compressor]=compressor_power_expressions
    print("objective function expressions set\n")
    # y ( Ax + Bx̄ + z), ||Ax + Bx̄ + z||^2
    L1 =[] # yt-1 ( Axt + Bx̄t + zt-1)
    L2 =[] # ρ/2 (||Axt + Bx̄t + zt||^2)
    Lλ= [] # λz
    Lz= [] # β/2 |z|^2
    sse=0.0
    for z in 1:l
        for i in part[1][z][:junction]
            push!(L1, JuMP.@NLexpression(m,(p[z][i]- p[l+1][i]+ Zp[z][i])*Yp[z][i] ))
            push!(L2, JuMP.@NLexpression(m,(p[z][i]- p[l+1][i]+ Zp[z][i])^2 ))
            push!(Lλ, JuMP.@NLexpression(m, Zp[z][i]*expr[z][:λp][i]))
            push!(Lz,  JuMP.@NLexpression(m,Zp[z][i]^2))
            sse+=(expr[z][:Zp][i])^2
        end
        for i in part[1][z][:pipe]
            push!(L1, JuMP.@NLexpression(m,(phi[z][i]- phi[l+1][i]+ Zphi[z][i])*Yphi[z][i] ))
            push!(L2, JuMP.@NLexpression(m,(phi[z][i]- phi[l+1][i]+ Zphi[z][i])^2 ))
            push!(Lλ, JuMP.@NLexpression(m, Zphi[z][i]*expr[z][:λphi][i]))
            push!(Lz,  JuMP.@NLexpression(m,Zphi[z][i]^2))
            sse+=(expr[z][:Zphi][i])^2
        end
        for i in part[1][z][:receipt]
            push!(L1, JuMP.@NLexpression(m,(supply[z][i]- supply[l+1][i]+ Zsupply[z][i])*Ysupply[z][i]))
            push!(L2, JuMP.@NLexpression(m,(supply[z][i]- supply[l+1][i]+ Zsupply[z][i])^2))
            push!(Lλ, JuMP.@NLexpression(m, Zsupply[z][i]*expr[z][:λs][i]))
            push!(Lz,  JuMP.@NLexpression(m,Zsupply[z][i]^2))
            sse+=(expr[z][:Zs][i])^2
        end
        for i in part[1][z][:delivery]
            push!(L1, JuMP.@NLexpression(m,(demand[z][i]-demand[l+1][i]+ Zdemand[z][i])*Ydemand[z][i]))
            push!(L2, JuMP.@NLexpression(m,(demand[z][i]-demand[l+1][i]+ Zdemand[z][i])^2))
            push!(Lλ, JuMP.@NLexpression(m, Zdemand[z][i]*expr[z][:λd][i]))
            push!(Lz,  JuMP.@NLexpression(m,Zdemand[z][i]^2))
            sse+=(expr[z][:Zd][i])^2
        end
        for i in part[1][z][:compressor]
            push!(L1, JuMP.@NLexpression(m,(alpha[z][i]- alpha[l+1][i]+ Zalpha[z][i])*Yalpha[z][i]))
            push!(L1, JuMP.@NLexpression(m,(flow[z][i]- flow[l+1][i]+ Zflow[z][i])*Yflow[z][i]) )
            push!(L2, JuMP.@NLexpression(m,(alpha[z][i]- alpha[l+1][i]+ Zalpha[z][i])^2 ))
            push!(L2, JuMP.@NLexpression(m,(flow[z][i]- flow[l+1][i]+ Zflow[z][i])^2 ))
            push!(Lλ, JuMP.@NLexpression(m, Zalpha[z][i]*expr[z][:λα][i]))
            push!(Lz,  JuMP.@NLexpression(m,Zalpha[z][i]^2))
            push!(Lλ, JuMP.@NLexpression(m, Zflow[z][i]*expr[z][:λf][i]))
            push!(Lz,  JuMP.@NLexpression(m,Zflow[z][i]^2))
            sse+=expr[z][:Zf][i]^2
            sse+=expr[z][:Zα][i]^2
        end
    end
    expr[1][:params][:L1] = L1
    expr[1][:params][:L2] = L2
    expr[1][:params][:Lλ] = Lλ
    expr[1][:params][:Lz] = Lz

    expr[1][:params][:tolerance]= 0.0
    expr[1][:params][:tolerancelimit]=0.0
    expr[1][:params][:outertolerance]=sqrt(sse)
    expr[1][:params][:objVal]=0.0
    print(sse, "  sum of squared errors expression values set\n")
    print(" global variables initial value\n")

    return (m, var, con, expr, part)

end
