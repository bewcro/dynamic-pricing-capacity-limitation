using Plots, DataFrames

lines_14 = CSV.read("Data/testcase/feeder15/lines.csv", DataFrame)
nodes_14 = CSV.read("Data/testcase/feeder15/nodes.csv",DataFrame)

grid14 = plot(grid = false, showaxis=false,xlims = (minimum(nodes_14.x)-1,maximum(nodes_14.x)+1), ylims = (minimum(nodes_14.y)-1,maximum(nodes_14.y)+1))
for i in 1:size(lines_14.node_f,1)
    plot!([nodes_14.x[nodes_14.index .== lines_14.node_f[i]][1],nodes_14.x[nodes_14.index .== lines_14.node_t[i]][1]],[nodes_14.y[nodes_14.index .== lines_14.node_f[i]][1],nodes_14.y[nodes_14.index .== lines_14.node_t[i]][1]],label=false,color="black")
end
scatter!(nodes_14.x,nodes_14.y, markersize = 15, label = false, color = "white")
for i in 1:size(nodes_14.index,1)
    annotate!(nodes_14.x[i],nodes_14.y[i], text(nodes_14.index[i],10))
end
scatter!([nodes_14.x[1]],[nodes_14.y[1]], markersize = 15, color = "red", label=false)
scatter!([nodes_14.x[1]],[nodes_14.y[1]+0.5], alpha=0,series_annotations = text("Substation", 10), color = "red", label=false)
savefig(grid14,"Figures/testcase14nodes.pdf")
display(plot!())


### 28 Grid ###
lines_28 = CSV.read("Data/testcase/feeder29/lines.csv",DataFrame)
nodes_28 = CSV.read("Data/testcase/feeder29/nodes.csv",DataFrame)

grid28 = plot(grid = false, showaxis=false,xlims = (minimum(nodes_28.x)-1,maximum(nodes_28.x)+1), ylims = (minimum(nodes_28.y)-1,maximum(nodes_28.y)+1))
for i in 1:size(lines_28.node_f,1)
    plot!([nodes_28.x[nodes_28.index .== lines_28.node_f[i]][1],nodes_28.x[nodes_28.index .== lines_28.node_t[i]][1]],[nodes_28.y[nodes_28.index .== lines_28.node_f[i]][1],nodes_28.y[nodes_28.index .== lines_28.node_t[i]][1]],label=false,color="black")
end
scatter!(nodes_28.x,nodes_28.y, markersize = 15, label = false, color = "white")
for i in 1:size(nodes_28.index,1)
    annotate!(nodes_28.x[i],nodes_28.y[i], text(nodes_28.index[i],10))
end
scatter!([nodes_28.x[1]],[nodes_28.y[1]], markersize = 15, color = "red", label=false)
scatter!([nodes_28.x[1]],[nodes_28.y[1]+0.6], alpha=0,series_annotations = text("Substation", 10), color = "red", label=false)
savefig(grid28,"Figures/testcase28nodes.pdf")
display(plot!())

### 56 Grid ###
lines_56 = CSV.read("Data/testcase/feeder57/lines.csv",DataFrame)
nodes_56 = CSV.read("Data/testcase/feeder57/nodes.csv",DataFrame)

grid56 = plot(grid = false, showaxis=false,xlims = (minimum(nodes_56.x)-1,maximum(nodes_56.x)+1), ylims = (minimum(nodes_56.y)-1,maximum(nodes_56.y)+1))
for i in 1:size(lines_56.node_f,1)
    plot!([nodes_56.x[nodes_56.index .== lines_56.node_f[i]][1],nodes_56.x[nodes_56.index .== lines_56.node_t[i]][1]],[nodes_56.y[nodes_56.index .== lines_56.node_f[i]][1],nodes_56.y[nodes_56.index .== lines_56.node_t[i]][1]],label=false,color="black",linewidth = 2)
end
scatter!(nodes_56.x,nodes_56.y, markersize = 15, label = false, color = "white")
for i in 1:size(nodes_56.index,1)
    annotate!(nodes_56.x[i],nodes_56.y[i], text(nodes_56.index[i],10))
end
scatter!([nodes_56.x[1]],[nodes_56.y[1]], markersize = 15, color = "red", label=false)
scatter!([nodes_56.x[1]],[nodes_56.y[1]+0.6], alpha=0,series_annotations = text("Substation", 10), color = "red", label=false)
savefig(grid56,"Figures/testcase56nodes.pdf")
display(plot!())



res = 200
price = CSV.read("Data/elspotprices.csv",DataFrame)
plot()
price = rand(1:24)
for i in b
    block = b * res/24*25/2
    plot()
end