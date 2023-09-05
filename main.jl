using LaTeXStrings, Statistics, ProgressBars, Gurobi, StatsPlots

################################
### LOADING SCRIPTS AND DATA ###
################################
println("Loading Scripts and Data...")

### Setting root and data path ###
root = dirname(@__FILE__)
data_path = joinpath(root, "Data/")

### Setting plotting theme ###
theme(:default)

### including scripts ###
include("model.jl")
include("data processing.jl")

### Pre-processing and importing data ###
data, node, line, ind_cost, tot_cost, ind_profile, tot_profile, residual, total_residual = dataprocessing_upperlevel(data_path)

println("------------------------")
println("--- Loading complete ---")
println("------------------------")

### Residual calculation ###
residual        = sum(eachcol(data["D"] .- data["PV"]))
total_residual  = sum(sum(eachcol(data["D"] .- data["PV"])))

### Running all models ###
#creating empty dictionaries to store results
full                = Dict()
full_uniform        = Dict()
full_proportional   = Dict()

#Creating lists with values we want to investigate for cap variation and tariff discount
betas   = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]
caps    = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]

for c in caps
    #cap and benchmark calculation for cap scenario c
    cap = cap_setting(total_residual,c,data["spot"])
    community_benchmark, individual_benchmark = benchmarks(ind_cost,ind_profile,tot_profile, cap)
    for b in betas
        println("-------------")
        println("BETA = ",b,", CAP = ",c)
        #Storing model results in the relevant dictionaries
        full[b,c]               = dynamic_pricing(data,node,line,individual_benchmark,b,cap,"full","none")
        full_uniform[b,c]       = dynamic_pricing(data,node,line,individual_benchmark,b,cap,"full","uniform")
        full_proportional[b,c]  = dynamic_pricing(data,node,line,individual_benchmark,b,cap,"full","proportional")
        #Checking complementarity of prosumer import and export
        if maximum(value.(full[b,c][:p⁺]).*value.(full[b,c][:p⁻])) > 0
            println("Complementarity doesn't hold for use case beta ",b,"and cap", c)
        elseif maximum(value.(full_uniform[b,c][:p⁺]).*value.(full_uniform[b,c][:p⁻])) > 0
            println("Complementarity doesn't hold for use case beta ",b,"and cap", c)
        elseif maximum(value.(full_proportional[b,c][:p⁺]).*value.(full_proportional[b,c][:p⁻])) > 0
            println("Complementarity doesn't hold for use case beta ",b,"and cap", c)
        end
        #Checking complementarity of battery charge and discharge
        if maximum(value.(full[b,c][:pᶜʰ]).*value.(full[b,c][:pᵈⁱˢ])) > 0
            println("Complementarity doesn't hold for use case beta ",b,"and cap", c)
        elseif maximum(value.(full_uniform[b,c][:pᶜʰ]).*value.(full_uniform[b,c][:pᵈⁱˢ])) > 0
            println("Complementarity doesn't hold for use case beta ",b,"and cap", c)
        elseif maximum(value.(full_proportional[b,c][:pᶜʰ]).*value.(full_proportional[b,c][:pᵈⁱˢ])) > 0
            println("Complementarity doesn't hold for use case beta ",b,"and cap", c)
        end
    end
end
println("----------------------------------------------")
println("Results Run and Complementarity Check Complete")
println("----------------------------------------------")

#########################################
### RESULTS PROCESSING AND STATISTICS ###
#########################################

### Objective value ###
objective = zeros((length(betas),length(caps)))
for b in eachindex(betas)
    for c in eachindex(caps)
        objective[b,c] = objective_value(full[betas[b],caps[c]])
    end
end

### Benefit ###
omega_plus                  = Dict()
omega_minus                 = Dict()
omega_plus_uniform          = Dict()
omega_minus_uniform         = Dict()
omega_plus_proportional     = Dict()
omega_minus_proportional    = Dict()
for b in betas
    for c in caps
        omega_plus[b,c]                 = Array(value.(full[b,c][:ω⁺]))
        omega_minus[b,c]                = Array(value.(full[b,c][:ω⁻]))
        omega_plus_uniform[b,c]         = Array(value.(full_uniform[b,c][:ω⁺]))
        omega_minus_uniform[b,c]        = Array(value.(full_uniform[b,c][:ω⁻]))
        omega_plus_proportional[b,c]    = Array(value.(full_proportional[b,c][:ω⁺]))
        omega_minus_proportional[b,c]   = Array(value.(full_proportional[b,c][:ω⁻]))
    end
end

### Prices ###
#empty dictionaries for prices
Price = Dict()
Price_Uniform = Dict()
Price_Proportional = Dict()

#Putting prices into a dataframe so we can call the summary statistics table
for b in betas
    for c in caps
        Price[b,c] = DataFrame(transpose(Array(value.(full[b,c][:x]))), ["Prosumer $i" for i in 1:14])
        Price_Uniform[b,c] = DataFrame(transpose(Array(value.(full_uniform[b,c][:x]))), ["Prosumer $i" for i in 1:14])
        Price_Proportional[b,c] = DataFrame(transpose(Array(value.(full_proportional[b,c][:x]))), ["Prosumer $i" for i in 1:14])
    end
end

################
### PLOTTING ###
################

### Illustration of Capacity Limitation Scenarios ###

cap_plot = plot(0:23,cap_setting(total_residual,1,data["spot"]), linetype=:steppost, label = "1", color =:red, xlim=(0,24), legend=:top, legend_columns=-1, xlabel = "Time of Day [h]", ylabel = "Capacity Limitation [kW]", yguidefontcolor=:red)
plot!([23,24],[cap_setting(total_residual,1,data["spot"])[end],cap_setting(total_residual,1,data["spot"])[end]], color =:red, label = false)
plot!(0:23, cap_setting(total_residual,0.5,data["spot"]), linetype=:steppost, linestyle=:dash, label = "0.5", color=:red)
plot!([23,24],[cap_setting(total_residual,0.5,data["spot"])[end],cap_setting(total_residual,0.5,data["spot"])[end]], linestyle=:dash, color =:red, label = false)
plot!(0:23,cap_setting(total_residual,0,data["spot"]), linetype=:steppost, linestyle=:dashdot, label = "0", color =:red)
plot!([23,24],[cap_setting(total_residual,0,data["spot"])[end],cap_setting(total_residual,0,data["spot"])[end]], linestyle=:dashdot, color =:red, label = false)
axis2 = twinx()
plot!(axis2,0:23,data["spot"],linetype=:steppost,yaxis="Spot Price [DKK/kWh]",color=:black,label=false, xlim = (0,24))
plot!(axis2,[23,24],[data["spot"][end],data["spot"][end]],linetype=:steppost,yaxis="Spot Price [DKK/kWh]",color=:black,label=false, xlim = (0,24))
savefig(cap_plot,"Figures/cap_vs_price.svg")

### Heatmap of Objective Values for different betas and caps ###
#making string for axis labels
beta_strings = []
cap_strings = []
for b in betas
    push!(beta_strings,string(b))
    push!(cap_strings,string(b))
end
#plotting heatmap of community cost
community_cost_plot = heatmap(cap_strings[5:11], beta_strings, objective[:,5:11],ylabel="Tariff Discount Factor",xlabel = "Capacity Limitation Variation Factor",seriescolor= cgrad(:RdYlGn_9,rev=true))
savefig(community_cost_plot,"Figures/community_cost.svg")


### Benefit Distribution Plots ###
distributionplot = groupedbar(1:14,[omega_minus[0.6,1].-omega_plus[0.6,1]  omega_minus_uniform[0.6,1].-omega_plus_uniform[0.6,1] omega_minus_proportional[0.6,1].-omega_plus_proportional[0.6,1]], xticks = 1:14, xlabel = "Community Member", ylabel = "Benefit compared to Benchmark [DKK]", label = ["No mechanism" "Equal" "Proportional"])
savefig(distributionplot,"Figures/bar_distribution.svg")

### Service Delivery Plots ###
residualplot = bar(0.5:23.5,residual,label = false,bar_width=1, xlims=(0,24), ylims = (-5,40), title = "Residual Demand", color =:blue, ylabel = "Import Power [kW]",left_margin = 6*Plots.mm, linecolor=:match)
plot!(0:23,cap_setting(total_residual,1,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,1,data["spot"])[end],cap_setting(total_residual,1,data["spot"])[end]],color=:black, label = false, linewidth = 2)

uncoordinated_plot = bar(0.5:23.5,tot_profile,label = false, bar_width = 1,xlims = (0,24),ylims = (-5,40), title = "Uncoordinated DR", color =:red, xlabel="Time of Day [h]", linecolor=:match)
plot!(0:23,cap_setting(total_residual,1,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,1,data["spot"])[end],cap_setting(total_residual,1,data["spot"])[end]],color=:black, label = false, linewidth = 2)

coordinated_plot = bar(0.5:23.5,value.(full[0.6,1][:pⁱᵐ]),label = false, title = "Coordinated DR", bar_width = 1,xlims = (0,24),ylims = (-5,40), color =:green, linecolor=:match)
plot!(0:23,cap_setting(total_residual,1,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,1,data["spot"])[end],cap_setting(total_residual,1,data["spot"])[end]],color=:black, label = false, linewidth = 2)

service_delivery_plot = plot(residualplot, uncoordinated_plot, coordinated_plot, layout = (1,3), size = (1200,300), bottom_margin = 6*Plots.mm)
savefig(service_delivery_plot,"Figures/service_delivery.svg")

### Cap Scenarios Plot ###
highplot = bar(0.5:23.5,value.(full[0.6,1][:pⁱᵐ]),label = false, title = "High Variation", bar_width = 1,xlims = (0,24),ylims = (0,11), color =:green, linecolor=:match, ylabel = "Import Power [kW]",left_margin = 6*Plots.mm)
plot!(0:23,cap_setting(total_residual,1,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,1,data["spot"])[end],cap_setting(total_residual,1,data["spot"])[end]],color=:black, label = false, linewidth = 2)

lowplot = bar(0.5:23.5,value.(full[0.6,0.5][:pⁱᵐ]),label = false, title = "Low Variation", bar_width = 1,xlims = (0,24),ylims = (0,11), color =:green, linecolor=:match, xlabel = "Time of Day [h]")
plot!(0:23,cap_setting(total_residual,0.5,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,0.5,data["spot"])[end],cap_setting(total_residual,0.5,data["spot"])[end]],color=:black, label = false, linewidth = 2)

constantplot = bar(0.5:23.5,value.(full[0.6,0][:pⁱᵐ]),label = false, title = "Constant", bar_width = 1,xlims = (0,24),ylims = (0,11), color =:green, linecolor=:match)
plot!(0:23,cap_setting(total_residual,0,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,0,data["spot"])[end],cap_setting(total_residual,0,data["spot"])[end]],color=:black, label = false, linewidth = 2)

cap_scenario_plot = plot(highplot, lowplot, constantplot, layout = (1,3), size = (1200,300), bottom_margin = 6*Plots.mm)
savefig(cap_scenario_plot,"Figures/cap_scenario.svg")