function partition1(partition, ref)
    # returns 3 things- the components in every paritition, the common ones and global number
    nw = ref[:nw][0]
    l = length(partition)
    junction = partition
    P = keys(nw[:pipe])
    C = keys(nw[:compressor])
    S= keys(nw[:dispatchable_receipt])
    D= keys(nw[:dispatchable_delivery])
    common = [Dict() for i in 1:l]
    component = [Dict() for i in 1:l]
    glob = Dict()
    glob[:junction]=[]
    glob[:pipe]=[]
    glob[:compressor]=[]
    glob[:delivery]=[]
    glob[:receipt]=[]

    for i in 1:l
        common[i][:pipe]=[]
        common[i][:junction]=[]
        common[i][:compressor] =[]
        common[i][:receipt]=[]
        common[i][:delivery]=[]
        component[i][:pipe]=[]
        component[i][:receipt]=[]
        component[i][:delivery] =[]
        component[i][:compressor] =[]
        component[i][:junction] = []
    end

    for z in 1:l
        for i in partition[z]
            p_fr = nw[:pipes_fr][i]
            p_to = nw[:pipes_to][i]
            c_fr=  nw[:compressors_fr][i]
            c_to=  nw[:compressors_to][i]
            for p in p_fr
                push!(component[z][:junction],nw[:pipe][p]["fr_junction"])
                push!(component[z][:junction],nw[:pipe][p]["to_junction"])
                push!(component[z][:pipe], p)
            end

            for p in p_to
                push!(component[z][:junction], nw[:pipe][p]["to_junction"])
                push!(component[z][:junction], nw[:pipe][p]["fr_junction"])
                push!(component[z][:pipe], p)
            end

            for c in c_fr
                push!(component[z][:junction],nw[:compressor][c]["fr_junction"])
                push!(component[z][:junction],nw[:compressor][c]["to_junction"])
                push!(component[z][:compressor], c)
            end

            for c in c_to
                push!(component[z][:junction], nw[:compressor][c]["to_junction"])
                push!(component[z][:junction], nw[:compressor][c]["fr_junction"])
                push!(component[z][:compressor], c)
            end
             for s in S
                 if ( nw[:dispatchable_receipt][s]["junction_id"] ==i)
                     push!(component[z][:receipt], s)
                 end
             end
             for d in D
                 if (nw[:dispatchable_delivery][d]["junction_id"] ==i)
                     push!(component[z][:delivery], d)
                 end
             end
        end
        component[z][:junction]= unique(component[z][:junction])
        component[z][:pipe]= unique(component[z][:pipe])
        component[z][:compressor]= unique(component[z][:compressor])
        component[z][:receipt]= unique(component[z][:receipt])
        component[z][:delivery]= unique(component[z][:delivery])
        common[z][:junction]=[i for i in component[z][:junction]if !(i in partition[z])]
        for i in common[z][:junction]
            for p in component[z][:pipe]
                if ((nw[:pipe][p]["fr_junction"] in partition[z]) &&  (nw[:pipe][p]["to_junction"] in partition[z]))
                    continue
                end
                push!(common[z][:pipe], p)
            end
            for c in component[z][:compressor]
                if ((nw[:compressor][c]["fr_junction"] in partition[z]) &&  (nw[:compressor][c]["to_junction"] in partition[z]))
                    continue
                end
                push!(common[z][:compressor], c)
            end
             for s in S
                 if ( nw[:dispatchable_receipt][s]["junction_id"] ==i)
                     push!(component[z][:receipt], s)
                     push!(common[z][:receipt],s)
                 end
             end
             for d in D
                 if ( nw[:dispatchable_delivery][d]["junction_id"] ==i)
                     push!(component[z][:delivery], d)
                     push!(common[z][:delivery],d)
                 end
             end
            common[z][:compressor]= unique(common[z][:compressor])
            common[z][:pipe]= unique(common[z][:pipe])
            component[z][:receipt] = unique(component[z][:receipt])
            common[z][:receipt] = unique(common[z][:receipt])
            component[z][:delivery] = unique(component[z][:delivery])
            common[z][:delivery]= unique(common[z][:delivery])

        end


    end
    for z in 1:l
        for h in 1: length(common[z][:junction])
            push!(glob[:junction], common[z][:junction][h])
        end
        for h in 1: length(common[z][:pipe])
            push!(glob[:pipe], common[z][:pipe][h])
        end
        for h in 1: length(common[z][:compressor])
            push!(glob[:compressor], common[z][:compressor][h])
        end
        for h in 1: length(common[z][:receipt])
            push!(glob[:receipt], common[z][:receipt][h])
        end
        for h in 1: length(common[z][:delivery])
            push!(glob[:delivery], common[z][:delivery][h])
        end
    end
    glob[:junction]= unique(glob[:junction])
    glob[:pipe]= unique(glob[:pipe])
    glob[:compressor]= unique(glob[:compressor])
    glob[:delivery]= unique(glob[:delivery])
    glob[:receipt]= unique(glob[:receipt])

return (component, common, glob)

end
