using CSV, DataFrames, Dates, Query, Plots, DataStructures, ProgressBars, LaTeXStrings, Statistics, ProgressBars, Gurobi, StatsPlots, Plots

root = dirname(@__FILE__)
grid_data = joinpath(root, "Data/European_LV_Test_Feeder_v2/European_LV_CSV/")

load_profiles_min = Dict()

D = zeros(55,24)
for i in 1:size(loads_data,1)
    load_profiles_min[i] = CSV.read(grid_data*"Load Profiles/Load_profile_$i.csv",DataFrame)

    for j in 1:size(load_profiles_min[i],1)
        if j == 1440
            D[i,1] += load_profiles_min[i].mult[j]/60
        else
            k = hour(Time(load_profiles_min[1].time[j]))
            D[i,k+1] += load_profiles_min[i].mult[j]/60
        end
    end
end

busses = CSV.read(grid_data*"Buscoords.csv",DataFrame,header=2)
line_data = CSV.read(grid_data*"Lines.csv",DataFrame,header=2)
line_code = CSV.read(grid_data*"LineCodes.csv",DataFrame,header=2)
loads_data = CSV.read(grid_data*"Loads.csv",DataFrame,header=3)


grid = scatter([busses." x"[1]],[busses." y"[1]],color="black", label = false, grid = false, showaxis = false)
for i in 1:length(line_data.Bus1)
        println([busses." x"[line_data.Bus1[i]],busses." x"[line_data.Bus2[i]]], [busses." y"[line_data.Bus1[i]],busses." y"[line_data.Bus2[i]]])
        plot!([busses." x"[line_data.Bus1[i]],busses." x"[line_data.Bus2[i]]], [busses." y"[line_data.Bus1[i]],busses." y"[line_data.Bus2[i]]],color = "black", label = false)
end
for b in loads_data.Bus
    if b % 1 == 0
    scatter!([busses." x"[b]],[busses." y"[b]], color = "red", markersize = 3, label = false)
    end
end
display(plot!())
savefig("Figures/grid case study.pdf")



node = Dict()
line = Dict()
for i in 1:size(line_data,1)
    ind = i
    node_f = line_data[i,:Bus1]
    node_t = line_data[i,:Bus2]
    print(line_data[i,:LineCode])
    print(line_data[i,:Length])
    r = line_data[i,:Length]/1000 * line_code.R1[line_code.Name .== line_data[i,:LineCode]]
    x = line_data[i,:Length]/1000 * line_code.X1[line_code.Name .== line_data[i,:LineCode]]
    f̅ = line_data[i,:s_max]
    add_line = Line(ind,node_f,node_t,r,x,f̅)
    line[add_line.ind] = add_line
end

line_code.R1[line_data.LineCode[1] .== line_code.Name] 