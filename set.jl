using JuMP
function setVariables(nw, l, part, β)
    m = JuMP.Model()
    var = [Dict() for i in 1:l+1]
    con = [Dict() for i in 1:l+1]
    expr = [Dict() for i in 1:l+1]
    #p = [Dict() for i in 1: l+1]
# This module defines local ( x, y, z) global variables ̄x and sets initial values of the global variables
# also defines all possible constraints and variables used in the program at some point.
    supply=Dict(); demand= Dict();phi=Dict();flow=Dict();alpha=Dict();p =Dict();
    Sp =  Dict();Sphi=Dict();Sflow=Dict();Salpha=Dict();
    Yp =  Dict();Ysupply=Dict(); Ydemand= Dict();Yphi=Dict();Yflow=Dict();Yalpha=Dict();
    Zp =  Dict();Zsupply=Dict(); Zdemand= Dict();Zphi=Dict();Zflow=Dict();Zalpha=Dict();
    for i in 1:l+1
            var[i][:compressor_power] = Dict()
            var[i][:compressor_powerE] = Dict()
            var[i][:net_nodal_injection] = Dict()
            var[i][:net_nodal_edge_out_flow] = Dict()
            var[i][:net_nodal_injectionE] = Dict()
            var[i][:net_nodal_edge_out_flowE] = Dict()
            con[i][:slack_pressure] = Dict()
            con[i][:nodal_balance] = Dict()
            con[i][:nodal_balanceE] = Dict()
            con[i][:pipe_physics] = Dict()
            con[i][:pipe_physicsE] = Dict()
            con[i][:compressor_flow_dir] = Dict()
            con[i][:compressor_boost] = Dict()
            con[i][:compressor_power] = Dict()
            con[i][:compressor_powerE1] = Dict()
            con[i][:compressor_flow_dirE] = Dict()
            con[i][:compressor_boostE] = Dict()
            con[i][:compressor_powerE] = Dict()
            con[i][:equalityP] = Dict()
            con[i][:equalityIn] = Dict()
            con[i][:equalityOut] = Dict()
            con[i][:equalityPhi] = Dict()
            con[i][:equalityD] = Dict()
            con[i][:equalityS] = Dict()
            con[i][:equalityFlow] = Dict()
            con[i][:equalityAlpha] = Dict()
            con[i][:equalityPower] = Dict()
            con[i][:equalityS]= Dict()
            con[i][:equalityD]= Dict()
            expr[i][:params] = Dict() # 1= fx elements, :fx, L1 and L2, 2 =somethingelse if need be

    end
    varname = ["X", "Y", "Z", "λ"]
    varname2 = ["p","phi","s","d","f","α"]
    varname3 = ["junction", "pipe","receipt","delivery","compressor", "compressor"]
    for i in 1: length(varname)
        s1 = varname[i]
        for j in 1: length(varname2)
            keycon =Symbol(string(s1, varname2[j]))
            for z in 1: l+1
                con[z][keycon]= Dict()
                expr[z][keycon]=Dict()
            end
        end
    end

    # defining local variables for each partition 'i'
    for z in 1: l
        p[z] =
            var[z][:pressure] = JuMP.@variable(
                m,
                [i in part[1][z][:junction]],
                lower_bound = nw[:junction][i]["p_min"],
                upper_bound = nw[:junction][i]["p_max"],
                base_name = "pressure"
            )
            print("Size Local Pressure variables in partition ",z, "\t", length(var[z][:pressure]),"\n")
        supply[z] =
            var[z][:fs] = JuMP.@variable(
                m,
                [i in part[1][z][:receipt]],
                lower_bound = nw[:receipt][i]["injection_min"],
                upper_bound = nw[:receipt][i]["injection_max"],
                base_name = "s"
            )
            print("SIZE Local supply variables in partition ",z, "\t", length(supply[z]),"\n")
        demand[z] =
            var[z][:fd] = JuMP.@variable(
                m,
                [i in part[1][z][:delivery]],
                lower_bound = nw[:delivery][i]["withdrawal_min"],
                upper_bound = nw[:delivery][i]["withdrawal_max"],
                base_name = "d"
            )
            print("SIZE Local demand variables in partition ",z, "\t", length(var[z][:fd]),"\n")
        phi[z] =
            var[z][:pipe_flux] = JuMP.@variable(
                m,
                [i in part[1][z][:pipe]],
                lower_bound = nw[:pipe][i]["flux_min"],
                upper_bound = nw[:pipe][i]["flux_max"],
                base_name = "phi"
            )
            print("SIZE Local flux variables in partition ",z, "\t", length(var[z][:pipe_flux]),"\n")
        flow[z] =
            var[z][:compressor_flow] = JuMP.@variable(
                m,
                [i in part[1][z][:compressor]],
                lower_bound = nw[:compressor][i]["flow_min"],
                upper_bound = nw[:compressor][i]["flow_max"],
                base_name = "flow"
            )

            print("SIZE Local flow variables in partition ",z, "\t", length(var[z][:compressor_flow]),"\n")
        alpha[z] =
            var[z][:compressor_ratio] = JuMP.@variable(
                m,
                [i in part[1][z][:compressor]],
                lower_bound = nw[:compressor][i]["c_ratio_min"],
                upper_bound = nw[:compressor][i]["c_ratio_max"],
                base_name = "alpha"
            )
            print("SIZE Local alpha variables in partition ",z, "\t", length(var[z][:compressor_ratio]),"\n")
            print("Local LoopDone \n")
    end
print("Local Done \n")
    # DEFINING GLOBAL VARIABLES
    p[l+1]= var[l+1][:pressure] = JuMP.@variable(
        m,
        [i in keys(nw[:junction])],
        lower_bound = nw[:junction][i]["p_min"],
        upper_bound = nw[:junction][i]["p_max"],
        base_name = "pressure"
    )
    print("Number of global pressure variables ","\t", length(var[l+1][:pressure]),"\n")

    for i in keys(nw[:junction])
        #if !(i in part[3][:junction])
            expr[l+1][:Xp][i] = 0.5*(nw[:junction][i]["p_min"]+ nw[:junction][i]["p_max"])
        #end
    end
    phi[l+1] =
        var[l+1][:pipe_flux] = JuMP.@variable(
            m,
            [i in keys(nw[:pipe])],
            lower_bound = nw[:pipe][i]["flux_min"],
            upper_bound = nw[:pipe][i]["flux_max"],
            base_name = "phi"
        )
        print("Number of global flux variables ", "\t", length(var[l+1][:pipe_flux]),"\n")
    for i in keys(nw[:pipe])
        #if !(i in part[3][:pipe])
            expr[l+1][:Xphi][i] = 0.5*(nw[:pipe][i]["flux_min"]+ nw[:pipe][i]["flux_max"])
        #end
    end
    flow[l+1] =
        var[l+1][:compressor_flow] = JuMP.@variable(
            m,
            [i in keys(nw[:compressor])],
            lower_bound = nw[:compressor][i]["flow_min"],
            upper_bound = nw[:compressor][i]["flow_max"],
            base_name = "flow"
        )
        print("Number of global flow variables ", "\t", length(var[l+1][:compressor_flow]),"\n")
    for i in keys(nw[:compressor])
        #if !(i in part[3][:compressor])
            expr[l+1][:Xf][i] = 0.5*(nw[:compressor][i]["flow_min"]+ nw[:compressor][i]["flow_max"])
        #end
    end
    alpha[l+1] =
        var[l+1][:compressor_ratio] = JuMP.@variable(
            m,
            [i in keys(nw[:compressor])],
            lower_bound = nw[:compressor][i]["c_ratio_min"],
            upper_bound = nw[:compressor][i]["c_ratio_max"],
            base_name = "alpha"
        )
        print("Number of global alpha variables ", "\t", length(var[l+1][:compressor_ratio]),"\n")
    for i in keys(nw[:compressor])
        #if !(i in part[3][:compressor])
            expr[l+1][:Xα][i] = 0.5*(nw[:compressor][i]["c_ratio_min"]+ nw[:compressor][i]["c_ratio_max"])
        #end
    end
    demand[l+1] =
        var[l+1][:fd] = JuMP.@variable(
            m,
            [i in keys(nw[:delivery])],
            lower_bound = nw[:delivery][i]["withdrawal_min"],
            upper_bound = nw[:delivery][i]["withdrawal_max"],
            base_name = "d"
        )
        print("Number of global demand variables ", "\t", length(var[l+1][:fd]),"\n")
    for i in keys(nw[:delivery])
            expr[l+1][:Xd][i] = 0.95*(nw[:delivery][i]["withdrawal_min"]+ nw[:delivery][i]["withdrawal_max"])
    end
    supply[l+1] =
        var[l+1][:fs] = JuMP.@variable(
            m,
            [i in keys(nw[:receipt])],
            lower_bound = nw[:receipt][i]["injection_min"],
            upper_bound = nw[:receipt][i]["injection_max"],
            base_name = "s"
        )
        print("Number of global supply variables ", "\t", length(var[l+1][:fs]),"\n")
    for i in keys(nw[:receipt])
            expr[l+1][:Xs][i] = 0.95*(nw[:receipt][i]["injection_min"]+ nw[:receipt][i]["injection_max"])
    end

 # creating Y and Z variables

    for z in 1: l
        Yp[z] =
            var[z][:pressureY] = JuMP.@variable(
                m,
                [i in part[1][z][:junction] ],
                lower_bound = -typemax(Float64),
                upper_bound = typemax(Float64),
                base_name = "pressureY"
            )
            print("Number of Y Pressure variables in partition ",z, "\t", length(var[z][:pressureY]),"\n")

        Zp[z] =
            var[z][:pressureZ] = JuMP.@variable(
                m,
                [i in part[1][z][:junction] ],
                lower_bound = -typemax(Float64),
                upper_bound = typemax(Float64),
                base_name = "pressureZ"
            )
            print("Number of Z Pressure variables in partition ",z, "\t", length(var[z][:pressureZ]),"\n")
        Ysupply[z] =
            var[z][:fsY] = JuMP.@variable(
                m,
                [i in part[1][z][:receipt]],
                lower_bound = -typemax(Float64),
                upper_bound = typemax(Float64),
                base_name = "sY"
        )
        print("Number of Y supply variables in partition ", "\t", length(var[z][:fsY]),"\n")
        Zsupply[z] =
            var[z][:fsZ] = JuMP.@variable(
                m,
                [i in part[1][z][:receipt]],
                lower_bound = -typemax(Float64),
                upper_bound = typemax(Float64),
                base_name = "sZ"
            )
            print("Number of Z supply variables in partiotoin", "\t", length(var[z][:fsZ]),"\n")
        Ydemand[z] =
            var[z][:fdY] = JuMP.@variable(
                m,
                [i in part[1][z][:delivery]],
                lower_bound = -typemax(Float64),
                upper_bound = typemax(Float64),
                base_name = "dY"
            )
            print("Number of Y demand variables in partition ",z, "\t", length(var[z][:fdY]),"\n")
        Zdemand[z] =
            var[z][:fdZ] = JuMP.@variable(
                m,
                [i in part[1][z][:delivery]],
                lower_bound = -typemax(Float64),
                upper_bound = typemax(Float64),
                base_name = "dZ"
            )
            print("Number of Z demand variables in partition ",z, "\t", length(var[z][:fdZ]),"\n")
        Yphi[z] =
            var[z][:pipe_fluxY] = JuMP.@variable(
                m,
                [i in part[1][z][:pipe]],
                lower_bound = -Inf,
                upper_bound = Inf,
                base_name = "phi-Y"
            )
            print("Number of Y flux variables in partition ",z, "\t", length(var[z][:pipe_fluxY]),"\n")
        Zphi[z] =
            var[z][:pipe_fluxZ] = JuMP.@variable(
                m,
                [i in part[1][z][:pipe]],# &&!(i in part[2][z][:pipe])],
                lower_bound = -Inf,
                upper_bound = Inf,
                base_name = "phi-Z"
            )
            print("Number of Z flux variables in partition ",z, "\t", length(var[z][:pipe_fluxZ]),"\n")
        Yflow[z] =
            var[z][:compressor_flowY] = JuMP.@variable(
                m,
                [i in part[1][z][:compressor]],# &&!(i in part[2][z][:compressor])],
                lower_bound = -Inf,
                upper_bound = Inf,
                base_name = "flow-Y"
            )
            print("Number of Y flow variables in partition ",z, "\t", length(var[z][:compressor_flowY]),"\n")
        Zflow[z] =
            var[z][:compressor_flowZ] = JuMP.@variable(
                m,
                [i in part[1][z][:compressor]],# &&!(i in part[2][z][:c ompressor])],
                lower_bound = -Inf,
                upper_bound = Inf,
                base_name = "flow-Z"
            )
            print("Number of Z flow variables in partition ",z, "\t", length(var[z][:compressor_flowZ]),"\n")
        Yalpha[z] =
            var[z][:compressor_ratioY] = JuMP.@variable(
                m,
                [i in part[1][z][:compressor]],# &&!(i in part[2][z][:compressor])],
                lower_bound = -Inf,
                upper_bound =Inf,
                base_name = "alpha-Y"
            )
            print("Number of Y alpha variables in partition ",z, "\t", length(var[z][:compressor_ratioY]),"\n")
        Zalpha[z] =
            var[z][:compressor_ratioZ] = JuMP.@variable(
                m,
                [i in part[1][z][:compressor]],# &&!(i in part[2][z][:compressor])],
                lower_bound = -Inf,
                upper_bound =Inf,
                base_name = "alpha-Z"
            )
            print("Number of Z alpha variables in partition ",z, "\t", length(var[z][:compressor_ratioZ]),"\n")
            c=0
            for p in 1: length(varname3)
                for i in part[1][z][Symbol(string(varname3[p]))]
                    c+=1
                    expr[z][Symbol(string(varname[2], varname2[p]))][i] =0.0
                    expr[z][Symbol(string(varname[3], varname2[p]))][i] =-1.0
                    expr[z][Symbol(string(varname[4], varname2[p]))][i] = β
                end
            end
                print("set", c," initial values of Y Z  and λ variables in expr \n")
    end
    for z in 1:l
        for (k,v) in expr[z]
            print(k,"\n")
            for (k2,v2) in expr[z][k]
                print(k2, "\t", v2, "\n")
            end
        end
    end


    arr=[p,supply,demand,phi,flow,alpha,Yp,Ysupply,Ydemand,Yphi,Yflow,Yalpha, Zp,Zsupply,Zdemand,Zphi,Zflow,Zalpha]

return m, var, con, expr, arr
end
