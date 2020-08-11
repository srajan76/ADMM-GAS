using JuMP
using Ipopt
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
# model = ( m, var, con, expr, part)
# m, var, con, expr, part === model[1], model[2], model[3], model[4] and model[5]
function L(model, ref, ρ,β, s)
        econ_weight = ref[:economic_weighting]
        nw = ref[:nw][0]
        l = length(model[5][1])
        # tried adding the original augmented Lagrangian of the whole problem, the λz + β/2|z|^2 is set to zero
        if s== "z"
                multiplier= [0,1,0]
        elseif s =="global"
                multiplier = [0,1,0]
        else
                multiplier= [1,1,0]
                # adding global ==local constraints only when optimizing local variables
                for z in 1:l
                        for j in model[5][2][z][:junction]
                                model[3][z][:equalityP][j]=@constraint(model[1], model[2][z][:pressure][j] ==model[2][l+1][:pressure][j])
                        end
                        for j in  model[5][2][z][:pipe]
                                 model[3][z][:equalityPhi][j] =  @constraint(model[1], model[2][z][:pipe_flux][j] ==model[2][l+1][:pipe_flux][j])
                        end
                        for j in  model[5][2][z][:compressor]
                                 model[3][z][:equalityAlpha][j] =  @constraint(model[1], model[2][z][:compressor_ratio][j] ==model[2][l+1][:compressor_ratio][j])
                                 model[3][z][:equalityFlow][j] =  @constraint(model[1], model[2][z][:compressor_flow][j] ==model[2][l+1][:compressor_flow][j])
                        end
                        for j in  model[5][3][:delivery]#nw[:delivery]
                                if j in model[5][2][z][:delivery]
                                        model[3][z][:equalityD][j] =  @constraint(model[1], model[2][z][:fd][j] ==model[2][l+1][:fd][j])
                                end
                        end
                         for j in  model[5][3][:receipt]#nw[:receipt]
                                 if j in model[5][2][z][:receipt]
                                         model[3][z][:equalityS][j] =  @constraint(model[1], model[2][z][:fs][j] ==model[2][l+1][:fs][j])
                                 end
                        end
                end
        end

        conname = ["X", "Y", "Z"]
        conname2= ["p","phi","s","d","f","α"]
        varname = ["pressure","pipe_flux","fs","fd","compressor_flow","compressor_ratio"]
        varname2= ["junction","pipe","receipt","delivery","compressor","compressor"]
        varnameInd= ["junction","pipe","dispatchable_receipt","dispatchable_delivery","compressor","compressor"]
        varname3 = ["","Y","Z"]
        keyvar = [Symbol(string(w,q)) for q in varname3 for w in varname]
        keycon = [Symbol(string(a,v)) for a in conname for v in conname2]
        keyind= repeat([Symbol(a) for a in varname2], length(conname))
        keyind2 = repeat([Symbol(a) for a in varnameInd], length(varname3))
        # I=[]
        # for p in 1:length(varname2)
        #         arr=[]
        #         if p==3 ||p==4
        #                 arr= [ i for i in keys(nw[Symbol(varnameInd[p])])]
        #                 push!(I,arr)
        #                 continue
        #         end
        #         arr = [ i for i in keys(nw[Symbol(varnameInd[p])]) if !( i in model[5][3][Symbol(varname2[p])]) ]
        #         push!(I,arr)
        # end





        # Adding constraints equating all the other variables other than the one being optimized being set
        # to their previous values
        for z in 1:l
                for p in 1:length(keycon)
                        if length(model[4][z][keycon[p]])==0
                                continue
                        end
                        for i in model[5][1][z][keyind[p]]
                                #if !(i  in model[5][2][z][keyind[p]] )
                                model[3][z][keycon[p]][i] = @constraint(
                                        model[1], model[2][z][keyvar[p]][i]== model[4][z][keycon[p]][i])
                                #end
                        end
                end
        end
        # If global variables are to be set to previous optimized value. doing it here
        for p in 1:length(varname2)

                if length(model[4][l+1][keycon[p]])==0
                        continue
                end

                for i in keys(nw[Symbol(varnameInd[p])])
                        model[3][l+1][keycon[p]][i] = @constraint(model[1],
                                model[2][l+1][keyvar[p]][i]== model[4][l+1][keycon[p]][i] )
                end
        end
        #print("opt started\n")
        JuMP.@NLobjective(
        model[1],
        Min,
        multiplier[1]* (econ_weight *sum(model[4][1][:params][:load][i] for i = 1:length(model[4][1][:params][:load])) +
        (1 - econ_weight) *sum(model[4][1][:params][:compressor][i] for i = 1:length(model[4][1][:params][:compressor]))) +
        multiplier[2]*(sum(model[4][1][:params][:L1][i] for i =1: length(model[4][1][:params][:L1]))+
        sum(ρ/2*model[4][1][:params][:L2][i] for i in 1: length(model[4][1][:params][:L2])))
        +multiplier[3]*(sum(model[4][1][:params][:Lλ][i] for i =1: length(model[4][1][:params][:Lλ]))
        +sum(β/2*model[4][1][:params][:Lz][i] for i =1: length(model[4][1][:params][:Lz])))
        )
        JuMP.set_optimizer(model[1], ipopt)
        JuMP.optimize!(model[1])
        #print("opt completed \n")
        if s =="local"
                model[4][1][:params][:objVal]=JuMP.objective_value(model[1])
        end
        # writing the solution in the expr dictionary
        for z in 1: l
                for p in 1:length(keycon)
                        if s=="z"
                                for i in model[5][1][z][keyind[p]]
                                        model[4][z][keycon[p]][i] = JuMP.value(model[2][z][keyvar[p]][i])
                                end
                        else
                                for i in model[5][1][z][keyind[p]]
                                        #if (( i in model[5][2][z][keyind[p]]) && string(keycon[p])[1]=='X')
                                        #continue
                                        #end
                                        model[4][z][keycon[p]][i] = JuMP.value(model[2][z][keyvar[p]][i])

                                end
                        end
                end
        end
        for p in 1:length(conname2)
                for i in keys(nw[Symbol(varnameInd[p])])
                        model[4][l+1][keycon[p]][i] = JuMP.value(model[2][l+1][keyvar[p]][i])
                end
        end
        #Printing global variables
        print("Printing variables After Optimization \n")
        print(" local\n ")
        for z in 1: l
                for p in 1:length(keycon)
                        for (i,v) in (model[4][z][keycon[p]])
                                print(keyvar[p], "\t", i,"\t",v, "\n")
                        end
                end
        end
        print(" global\n ")
        for p in 1:length(varname)
                for (i,v) in (model[4][l+1][keycon[p]])
                        print(keyvar[p], "\t", i,"\t",v, "\n")
                end
        end

        #printVal(model)
        # deleting constraints by var instead of removing the constraint dictionary
        for z in 1:l+1
                for p in 1:length(keycon)
                        if length(model[3][z][keycon[p]])==0
                                continue
                        end
                        for i in keys(model[3][z][keycon[p]])
                                #print( model[5][z][keyind[p]])
                                #print("I", i,"Z",z,"P",keycon[p], "\n")
                                delete(model[1], model[3][z][keycon[p]][i])
                        end
                end
        end
        conX= ["equalityP","equalityPhi","equalityAlpha","equalityFlow", "equalityD","equalityS"]
        for z in 1:l
                for j in conX
                        if length(model[3][z][Symbol(j)])==0
                                continue
                        end
                        print(keys(model[3][z][Symbol(j)]))
                        for i in keys(model[3][z][Symbol(j)])
                                delete(model[1], model[3][z][Symbol(j)][i])
                        end
                end
        end
        for z in 1:l
                for j in conX
                        if( length(model[3][z][Symbol(string(j))])>0)
                                for (k,v) in model[3][z][Symbol(string(j))]
                                        pop!(model[3][z][Symbol(string(j))],k)
                                end
                        end
                end
        end


        for z in 1:l+1
                for i in conname
                        for j in conname2
                                if( length(model[3][z][Symbol(string(i,j))])>0)
                                        for (k,v) in model[3][z][Symbol(string(i,j))]
                                                pop!(model[3][z][Symbol(string(i,j))],k)
                                        end
                                end
                        end
                end
        end
        #print(" OPT done")
return model
end

function optOver!(A, model)
        #print("enter optover", "\n")
        conname = ["X", "Y", "Z"]
        conname2 = ["p","phi","s","d","f","α"]
        varname= ["junction","pipe","receipt","delivery","compressor","compressor"]

        l = length(model[5][1])
        if A=="Global"
                i = conname[1]
                for j in conname2
                        model[4][l+1][Symbol(string(i,j))]= Dict()
                end
                # also remove the xlocal which are shared
                # for z in 1:l
                #         for j in 1: length(conname2)
                #                 if j==3 || j==4
                #                         continue
                #                 end
                #                 for (k,v) in model[4][z][Symbol(string(i,conname2[j]))]
                #                         if k in (model[5][2][z][Symbol(varname[j])])
                #                                 pop!(model[4][z][Symbol(string(i,conname2[j]))],k)
                #                         end
                #                 end
                #         end
                # end
        elseif A =="X"
                i=A
                # for j in 1: length(conname2)
                #         for z in 1:l
                #                 model[4][z][Symbol(string(i,conname2[j]))]= Dict()
                #         end
                #         if j ==3||j==4
                #                 continue
                #         end
                #         for (k,v) in model[4][l+1][Symbol(string(i,conname2[j]))]
                #                 if k in model[5][3][Symbol(varname[j])]
                #                         pop!(model[4][l+1][Symbol(string(i,conname2[j]))],k)
                #                 end
                #         end
                # end
                #i = conname[1]
                for z in 1:l
                        for j in conname2
                                model[4][z][Symbol(string(i,j))]= Dict()
                        end
                end

        else
                i =A
                for j in conname2
                        for z in 1:l
                                model[4][z][Symbol(string(i,j))]= Dict()
                        end
                end
        end
# if X then remove also xglobal whch are shared
# if Z it is fine
end


function updateY!(model, ref,ρ,k)
        nw = ref[:nw][0]
        l = length(model[5][1])
        varname = ["pressure","pipe_flux","fs","fd","compressor_flow","compressor_ratio"]
        varname2= ["junction","pipe","receipt","delivery","compressor","compressor"]
        varname3 = ["","Y","Z"]
        conname= ["X","Y","Z"]
        conname2 = ["p","phi","s","d","f","α"]
        varnameInd= ["junction","pipe","dispatchable_receipt","dispatchable_delivery","compressor","compressor"]
        keyvar = [Symbol(string(w,q)) for q in varname3 for w in varname]

        keyind= [Symbol(a) for a in varname2]
        keycon =[[Symbol(string(a,v)) for v in conname2] for a in conname]

        y = [Dict() for i in 1:l]
        sz=0
        for z in 1:l
                for j in conname2
                        y[z][Symbol(string("Y", j))]= model[4][z][Symbol(string("Y", j))]
                        sz+= length(model[4][z][Symbol(string("Y", j))])
                end
        end
        tol=0.0
        tollim= 2*sqrt(sz)/(k*ρ)
        #tollim= sqrt(sz)/1000
        #tollim= 0.05
        # update y(t+1) = y(t) + ρ( Axt + Bx̄t + zt) and calculate ||Axt + Bx̄t + zt||^2
        for z in 1:l
                for j in 1: length(conname2)
                        idx=[Symbol(string(c, conname2[j])) for c in conname]
                        for k in  model[5][1][z][keyind[j]]
                                tmp = model[4][z][idx[1]][k]- model[4][l+1][idx[1]][k]+ model[4][z][idx[3]][k]
                                 tol+= tmp^2
                                #tol+= (tmp- model[4][z][idx[3]][k])^2
                                model[4][z][idx[2]][k] = y[z][idx[2]][k] + ρ*(tmp)
                        end
                end
        end
        model[4][1][:params][:tolerance]=sqrt(tol)
        model[4][1][:params][:tolerancelimit]= tollim
        # Printing global variables
        print("Printing variables After Optimization \n")
        print(" local\n ")
        keycon = [Symbol(string(a,v)) for a in conname for v in conname2]
        for z in 1: l
                for p in 1:length(keycon)
                        for (i,v) in (model[4][z][keycon[p]])
                                print(keyvar[p], "\t", i,"\t",v, "\n")
                        end
                end
        end
        print(" global\n ")
        for p in 1:length(varname)
                for (i,v) in (model[4][l+1][keycon[p]])
                        print(keyvar[p], "\t", i,"\t",v, "\n")
                end
        end


end
function updateλ!(model,β, λ_h, λ_l)
        conname2= ["p","phi","s","d","f","α"]
        varname2= ["junction","pipe","receipt","delivery","compressor","compressor"]
        l = length(model[5][1])
        lamb=[Dict() for z in 1:l]
        for z in 1:l
                for j in conname2
                        lamb[z][Symbol(string("λ",j))]= model[4][z][Symbol(string("λ",j))]
                end
        end

        dot =0.0
        counter =0
        sse=0.0
        for z in 1:l
                for j in conname2
                        counter+= length(lamb[z][Symbol(string("λ",j))])
                        if length(lamb[z][Symbol(string("λ",j))])>0
                                for(k,v) in lamb[z][Symbol(string("λ",j))]
                                        model[4][z][Symbol(string("λ",j))][k]= v+β*model[4][z][Symbol(string("Z",j))][k]
                                        dot+= (λ_h - λ_l)*model[4][z][Symbol(string("λ",j))][k]
                                        sse+=model[4][z][Symbol(string("Z",j))][k]^2
                                end
                        end
                end
        end
        for z in 1:l
                for j in conname2
                        for(k,v) in lamb[z][Symbol(string("λ",j))]
                                model[4][z][Symbol(string("λ",j))][k]*=dot/(counter*((λ_h-λ_l)^2))
                        end
                end
        end



        model[4][1][:params][:outertolerance]=sqrt(sse)


        #printing Lambda
        # print("final λ \n")
        # for z in 1:l
        #         for j in conname2
        #                 idx = Symbol(string("λ",j))
        #                 for (k,v) in model[4][z][idx]
        #                         print(idx, "\t", k, "\t", v, "\n")
        #                 end
        #         end
        # end

end

function reinitialize!( model, β,ref)
        conname2= ["p","phi","s","d","f","α"]
        varname2= ["junction","pipe","receipt","delivery","compressor","compressor"]
        conname= ["Y","λ","Z","X"]
        nw = ref[:nw][0]
        l = length(model[5][1])
        #  setting  λ+ βz + y =0 ⟹ z= -1/β, y = 1- λ
        for z in 1:l
                for j in 1: length(conname2)

                        for (k,v) in model[4][z][Symbol(string(conname[3],conname2[j]))]
                                model[4][z][Symbol(string(conname[3],conname2[j]))][k] = -1/β
                        end
                        for (k,v) in model[4][z][Symbol(string(conname[1],conname2[j]))]
                                model[4][z][Symbol(string(conname[1],conname2[j]))][k] = 1-model[4][z][Symbol(string(conname[2],conname2[j]))][k]
                        end
                end
        end
        # tried reinitializingg  ̄x.. not cnverging..
        # for i in keys(nw[:junction])
        #         model[4][l+1][:Xp][i] = 0.5*(nw[:junction][i]["p_min"]+ nw[:junction][i]["p_max"])
        # end
        # for i in keys(nw[:pipe])
        #         model[4][l+1][:Xphi][i] = 0.5*(nw[:pipe][i]["flux_min"]+ nw[:pipe][i]["flux_max"])
        # end
        # for i in keys(nw[:compressor])
        #     model[4][l+1][:Xf][i] = 0.5*(nw[:compressor][i]["flow_min"]+ nw[:compressor][i]["flow_max"])
        # end
        # for i in keys(nw[:compressor])
        #        model[4][l+1][:Xα][i] = 0.5*(nw[:compressor][i]["c_ratio_min"]+ nw[:compressor][i]["c_ratio_max"])
        # end
        # for i in keys(nw[:delivery])
        #        model[4][l+1][:Xd][i] = 0.5*(nw[:delivery][i]["withdrawal_min"]+ nw[:delivery][i]["withdrawal_max"])
        # end
        # for i in keys(nw[:receipt])
        #       model[4][l+1][:Xs][i] = 0.5*(nw[:receipt][i]["injection_min"]+ nw[:receipt][i]["injection_max"])
        # end
end

function printVal(model)
        l = length(model[5][1])
        for z in 1:l+1
                for (k,v) in model[4][z]
                        if length(model[4][z][k])>0
                                for (k2,v2) in model[4][z][k]
                                        print(k2, "\t", v2, "\n")
                                end
                        end
                end
        end
return
end
