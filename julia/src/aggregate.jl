############################################################################################################
################################################################################################ Julia logs
############################################################################################################

# draw a png plot of the a-norm and 2-norm errors, and then save these files locally in the same folder
function drawJuliaPlot(path::ASCIIString, treeIndex::ASCIIString; only2=false)
    names,dataA,data2 = aggJuliaRes(path, treeIndex);

    if names == -1
        return
    end

    #first plot julia default
    for i in 1:length(names)
        if search(names[i], "default").stop != -1
            plot(log(data2[i]), label=names[i], linewidth=2.5)
        end
    end

    for p2 in 6:15
        for i in 1:length(names)
            if search(names[i], string(2^p2)).stop != -1
                plot(log(data2[i]), label=names[i], linewidth=2.5)
            end
        end
    end
    legend(fontsize=7)
    ylabel("Relative 2-norm")
    xlabel("Iterations")

    savefig(path * treeIndex * "_2norm_julia.png")

    PyPlot.clf()

    if only2
        return
    end

    #first plot julia default
    for i in 1:length(names)
        if search(names[i], "default").stop != -1
            plot(log(dataA[i]), label=names[i], linewidth=2.5)
        end
    end

    for p2 in 6:15
        for i in 1:length(names)
            if search(names[i], string(2^p2)).stop != -1
                plot(log(dataA[i]), label=names[i], linewidth=2.5)
            end
        end
    end
    legend(fontsize=7)
    ylabel("Relative A-norm")
    xlabel("Iterations")

    savefig(path * treeIndex * "_Anorm_julia.png")

    PyPlot.clf()
end

# aggregates the results in julia format from the folder in 'path' for the tree 'treeIndex'
function aggJuliaRes(path::ASCIIString, treeIndex::ASCIIString)

    names = Array{ASCIIString,1}(0)
    dataA = Array{Array{Float64,1},1}(0)
    data2 = Array{Array{Float64,1},1}(0)

    haveData = false

    allLogs = readdir(path)
    for log in allLogs
        if search(log, "julia").stop != -1 && search(log, treeIndex).stop != -1 && search(log, "png").stop == -1
            curData = parseData(path * log)

            normA = Array{Float64,1}(0)
            norm2 = Array{Float64,1}(0)
            for i in 1:length(curData)
                push!(normA, curData[i][1])
                push!(norm2, curData[i][2])
            end

            name = "Julia_default"
            for i in 6:15
                if search(log, string(2^i)).stop != -1
                    name = "Julia_" * string(2^i)
                    break
                end
            end

            push!(names, name)
            push!(dataA, normA)
            push!(data2, norm2)

            haveData = true
        end
    end

    if !haveData
        return -1,-1,-1
    end

    return names, dataA, data2

end

# parse a julia log file, and return the anorm and 2norm errors (in this order)
function parseData(fn::ASCIIString; only2 = false)
    f = open(fn)
    lines = readlines(f)
    close(f)
    
    results = Array{Array{Float64,1},1}(0)

    nr = parse(Int64, split(lines[2], ' ')[1])
    for i in 3:(3 + nr - 1)
        wholeLn = split(lines[i], [' ', '=', '\n', '\t'])
        ln = []
        for j in wholeLn
            if j != ""
                push!(ln, j)
            end
        end
        
        thisLine = []
        for j in 4:2:length(ln)
            push!(thisLine, parse(Float64, ln[j]))
        end
        
        push!(results, thisLine)
    end

    # return only the two norm
    if only2
        results2 = []
        for i in 1:length(results)
            push!(results2, results[i][2])
        end
        results = copy(results2)
    end
    
    return results
    
end

# parse a julia log file, and return the anorm and 2norm errors
function parsePythonData(fn::ASCIIString)
    f = open(fn)
    lines = readlines(f)
    close(f)
    
    results = Array{Array{Float64,1},1}(0)

    nr = length(lines) - 1
    for i in 2:(nr + 1)
        wholeLn = split(lines[i], [','])

        ln = []
        for j in wholeLn
            if j != ""
                push!(ln, parse(Float64, j))
            end
        end
        push!(results, ln)
    end
    
    return results
    
end

# parse a julia log file, and return the anorm and 2norm errors
function parseMatlabData(fn::ASCIIString)
    f = open(fn)
    lines = readlines(f)
    close(f)
    
    results = Array{Array{Float64,1},1}(0)

    nr = length(lines)
    for i in 7:nr
        wholeLn = split(lines[i], [' ', '=', ',', '\n'])
        ln = []
        for j in wholeLn
            if j != ""
                push!(ln, j)
            end
        end
        
        thisLine = []
        push!(thisLine, parse(Float64, ln[6]))
        push!(thisLine, parse(Float64, ln[4]))
        
        push!(results, thisLine)
    end
    
    return results
    
end


# parse a julia log file, and return the anorm and 2norm errors
function parseCPPData(fn::ASCIIString)
    f = open(fn)
    lines = readlines(f)
    close(f)
    
    results = Array{Float64,1}(0)

    nr = length(lines)
    for i in 2:nr-1
        wholeLn = split(lines[i], [' ', '=', ',', '\n', '\t'])
        ln = []
        for j in wholeLn
            if j != ""
                push!(ln, j)
            end
        end
        
        push!(results, parse(Float64, ln[end]))
    end
    
    return results
    
end
