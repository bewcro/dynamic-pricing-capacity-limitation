using LaTeXStrings, Statistics, ProgressBars, Gurobi, StatsPlots, Plots

#############################################
### LOADING SCRIPTS AND DATA MANIPULATION ###
#############################################
println("Loading Scripts and Data...")

### Setting root and data path ###
root = dirname(@__FILE__)
data_path = joinpath(root, "Data/")

### Setting plotting defaults ###
default(fontfamily = "Computer Modern", dpi = 400, size = (800, 600))

### including scripts ###
include("model.jl")
include("data processing.jl")

### Pre-processing and importing data ###
N_consumer = 14 #number of prosumers
grid = "feeder15" #grid to use
data, node, line, ind_cost, tot_cost, ind_profile, tot_profile = data_processing(data_path,grid,N_consumer)

println("Data loaded")


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
        full[b,c]               = dynamic_pricing(data,node,line,individual_benchmark,b,cap,"none",true)
        full_uniform[b,c]       = dynamic_pricing(data,node,line,individual_benchmark,b,cap,"uniform",true)
        full_proportional[b,c]  = dynamic_pricing(data,node,line,individual_benchmark,b,cap,"proportional",true)
    end
end

println("               Simulations Run                ")
println("----------------------------------------------")

#Checking complementarity of prosumer import and export
for c in caps
    for b in betas
        if maximum(value.(full[b,c][:p⁺]).*value.(full[b,c][:p⁻])) > 1e-6
            println("Power complementarity doesn't hold for no distribution use case beta ",b," and cap", c)
        end
        if maximum(value.(full_uniform[b,c][:p⁺]).*value.(full_uniform[b,c][:p⁻])) > 1e-6
            println("Power complementarity doesn't hold for uniform distribution use case beta ",b," and cap", c)
            println(maximum(value.(full_uniform[b,c][:p⁺]).*value.(full_uniform[b,c][:p⁻])))
        end
        if maximum(value.(full_proportional[b,c][:p⁺]).*value.(full_proportional[b,c][:p⁻])) > 1e-6
            println("Power complementarity doesn't hold for proportional distribution use case beta ",b," and cap", c)
        end
        #Checking complementarity of battery charge and discharge
        if maximum(value.(full[b,c][:pᶜʰ]).*value.(full[b,c][:pᵈⁱˢ])) > 1e-6
            println("Battery complementarity doesn't hold for no distribution use case beta ",b," and cap", c)
        end
        if maximum(value.(full_uniform[b,c][:pᶜʰ]).*value.(full_uniform[b,c][:pᵈⁱˢ])) > 1e-6
            println("Battery complementarity doesn't hold for uniform distribution use case beta ",b,"and cap", c)
        end
        if maximum(value.(full_proportional[b,c][:pᶜʰ]).*value.(full_proportional[b,c][:pᵈⁱˢ])) > 1e-6
            println("Battery complementarity doesn't hold for proportional distribution use case beta ",b,"and cap", c)
        end
    end
end

println("        Complementarity check complete        ")
println("----------------------------------------------")

#Checking if Demand shedding occurs

for c in caps
    for b in betas
        if maximum(value.(full[b,c][:dˢʰᵉᵈ])) >= 1e-6
            println("Demand shedding occurs in case study beta = $b, cap = $c")
        end
    end
end

println("        Demand shedding check complete        ")
println("----------------------------------------------")

#=
k=0
for c in caps
    for b in betas
        for t in 1:24
            for l in 1:14
                if value.(full[b,c][:f_p][l,t]/data["Sᵇᵃˢᵉ"])^2 + value.(full[b,c][:f_q][l,t]/data["Sᵇᵃˢᵉ"])^2 >= (line[l].f̅)^2
                    println("Case beta = $b and cap = $c with no distribution: Line $i capacity exceeded in hour $t by", value.(full[b,c][:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(full[b,c][:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 - (line[l].f̅)^2)
                    k += 1
                end
                if value.(full[b,c][:f_p][l,t]/data["Sᵇᵃˢᵉ"])^2 + value.(full[b,c][:f_q][l,t]/data["Sᵇᵃˢᵉ"])^2 >= (line[l].f̅)^2
                    println("Case beta = $b and cap = $c with no distribution: Line $i capacity exceeded in hour $t by", value.(full[b,c][:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(full[b,c][:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 - (line[l].f̅)^2)
                    k += 1
                end
            end
        end
    end
end
if k == 0
    println("No line capacities exceeded in any cases")
end

println("----------------------------------------------")
println("Ex-post line capacity check complete")
println("----------------------------------------------")
=#

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

describe(Price[0.6,1])
describe(Price_Uniform[0.6,1])
describe(Price_Proportional[0.6,1])

################
### PLOTTING ###
################
default(fontfamily = "Computer Modern", dpi = 600)

### Illustration of Capacity Limitation Scenarios ###
cap_plot = Plots.plot(0:23,cap_setting(total_residual,1,data["spot"]), linetype=:steppost, label = "1", color =:red, xlim=(0,24), legend=:top, legend_columns=-1, xlabel = "Time of Day [h]", ylabel = "Capacity Limitation [kW]", yguidefontcolor=:red)
plot!([23,24],[cap_setting(total_residual,1,data["spot"])[end],cap_setting(total_residual,1,data["spot"])[end]], color =:red, label = false)
plot!(0:23, cap_setting(total_residual,0.5,data["spot"]), linetype=:steppost, linestyle=:dash, label = "0.5", color=:red)
plot!([23,24],[cap_setting(total_residual,0.5,data["spot"])[end],cap_setting(total_residual,0.5,data["spot"])[end]], linestyle=:dash, color =:red, label = false)
plot!(0:23,cap_setting(total_residual,0,data["spot"]), linetype=:steppost, linestyle=:dashdot, label = "0", color =:red)
plot!([23,24],[cap_setting(total_residual,0,data["spot"])[end],cap_setting(total_residual,0,data["spot"])[end]], linestyle=:dashdot, color =:red, label = false)
axis2 = twinx()
plot!(axis2,0:23,data["spot"],linetype=:steppost,yaxis="Spot Price [DKK/kWh]",color=:black,label=false, xlim = (0,24))
plot!(axis2,[23,24],[data["spot"][end],data["spot"][end]],linetype=:steppost,yaxis="Spot Price [DKK/kWh]",color=:black,label=false, xlim = (0,24))
Plots.savefig(cap_plot,"Figures/cap_vs_price.pdf")

### Heatmap of Objective Values for different betas and caps ###
#making string for axis labels
beta_strings = []
cap_strings = []
for b in betas
    push!(beta_strings,string(b))
    push!(cap_strings,string(b))
end

community_cost_plot = Plots.heatmap(beta_strings, cap_strings, objective,xlabel="Capacity Limitation Variation Factor",ylabel = "Tariff Discount Factor",seriescolor= cgrad(:RdYlGn_9,rev=true), ylabelfontsize = 14, xlabelfontsize = 14, ytickfontsize=14,xtickfontsize=14)
Plots.savefig(community_cost_plot,"Figures/community_cost.pdf")

### Pricing Heatmaps ###
none = heatmap(1:24,1:14,transpose(Matrix(Price[0.5,0.5])),seriescolor=:heat, xlabel = "Time-of-day [h]", yticks = 1:14, ylabel = "Community member",left_margin= 10*Plots.mm, title = "No Distribution")
uniform = heatmap(1:24,1:14,transpose(Matrix(Price_Uniform[0.5,0.5])),seriescolor=:heat, xlabel = "Time-of-day [h]", yticks = false,left_margin= -10*Plots.mm, title = "Equal Distribution" #=, ylabel = "Community member"=#)
proportional = heatmap(1:24,1:14,transpose(Matrix(Price_Proportional[0.5,0.5])),seriescolor=:heat, xlabel = "Time-of-day [h]", yticks = false,left_margin= -10*Plots.mm, title = "Proportional Distribution"  #=, ylabel = "Community member"=#)
price_heatmaps = plot(none,uniform,proportional,layout = (1,3), size = (1600,400), ylabelfontsize=16,ytickfontsize = 14, xlabelfontsize=16,xtickfontsize=14,bottom_margin=10*Plots.mm,titlefontsize = 16, top_margin=3*Plots.mm)

savefig(price_heatmaps,"Figures/heatmap_prices.pdf")

### Benefit Distribution Plots ###
distributionplot = groupedbar(1:14,[omega_minus[0.5,0.5].-omega_plus[0.5,0.5] omega_minus_proportional[0.5,0.5].-omega_plus_proportional[0.5,0.5] omega_minus_uniform[0.5,0.5].-omega_plus_uniform[0.5,0.5]], xticks = 1:14, xlabel = "Community Member", ylabel = "Benefit compared to Benchmark [DKK]", label = ["No mechanism" "Proportional" "Equal"], legend = :outertop, legend_columns = 3, seriescolor = [:grey :purple :orange],xtickfontsize=12,ytickfontsize=12,xlabelfontsize=14,ylabelfontsize=14,legendfontsize=12)
Plots.savefig(distributionplot,"Figures/bar_distribution.pdf")

### Service Delivery Plots ###
residualplot = bar(0.5:23.5,residual,label = false,bar_width=1, xlims=(0,24), title = "Residual Demand", color =:blue, ylabel = "Import Power [kW]",left_margin = 8Plots.mm, linecolor=:match,ylabelfontsize=14,xlabelfontsize=14,ytickfontsize=14,xtickfontsize=14)
plot!(0:23,cap_setting(total_residual,0.5,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,0.5,data["spot"])[end],cap_setting(total_residual,0.5,data["spot"])[end]],color=:black, label = false, linewidth = 2)

uncoordinated_plot = bar(0.5:23.5,tot_profile,label = false, bar_width = 1,xlims = (0,24), title = "Uncoordinated DR", color =:red, xlabel="Time of Day [h]", linecolor=:match,ylabelfontsize=14,xlabelfontsize=14,ytickfontsize=14,xtickfontsize=14)
plot!(0:23,cap_setting(total_residual,0.5,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,0.5,data["spot"])[end],cap_setting(total_residual,0.5,data["spot"])[end]],color=:black, label = false, linewidth = 2)

coordinated_plot = bar(0.5:23.5,value.(full[0.5,0.5][:pⁱᵐ])-value.(full[0.5,1][:pᵉˣ]),label = false, title = "Coordinated DR", bar_width = 1,xlims = (0,24), color =:green, linecolor=:match,ylabelfontsize=14,xlabelfontsize=14,ytickfontsize=14,xtickfontsize=14)
plot!(0:23,cap_setting(total_residual,0.5,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,0.5,data["spot"])[end],cap_setting(total_residual,0.5,data["spot"])[end]],color=:black, label = false, linewidth = 2)

service_delivery_plot = Plots.plot(residualplot, uncoordinated_plot, coordinated_plot, layout = (1,3), size = (1200,300), bottom_margin = 10*Plots.mm, top_margin = 3*Plots.mm)
Plots.savefig(service_delivery_plot,"Figures/service_delivery.pdf")

### Cap Scenarios Plot ###
highplot = bar(0.5:23.5,value.(full[0.5,1][:pⁱᵐ]),label = false, title = "High Variation", bar_width = 1,xlims = (0,24),ylims = (0,3), color =:green, linecolor=:match, ylabel = "Import Power [kW]",left_margin = 8Plots.mm,ylabelfontsize=14,xlabelfontsize=14,ytickfontsize=12,xtickfontsize=12)
plot!(0:23,cap_setting(total_residual,1,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,1,data["spot"])[end],cap_setting(total_residual,1,data["spot"])[end]],color=:black, label = false, linewidth = 2)

lowplot = bar(0.5:23.5,value.(full[0.5,0.5][:pⁱᵐ]),label = false, title = "Low Variation", bar_width = 1,xlims = (0,24),ylims = (0,3), color =:green, linecolor=:match, xlabel = "Time of Day [h]",ylabelfontsize=14,xlabelfontsize=14,ytickfontsize=12,xtickfontsize=12)
plot!(0:23,cap_setting(total_residual,0.5,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,0.5,data["spot"])[end],cap_setting(total_residual,0.5,data["spot"])[end]],color=:black, label = false, linewidth = 2)

constantplot = bar(0.5:23.5,value.(full[0.5,0][:pⁱᵐ]),label = false, title = "Constant", bar_width = 1,xlims = (0,24),ylims = (0,3), color =:green, linecolor=:match,ylabelfontsize=14,xlabelfontsize=14,ytickfontsize=12,xtickfontsize=12)
plot!(0:23,cap_setting(total_residual,0,data["spot"]), linetype=:steppost, color=:black, xlims=(0,24),label=false, linewidth = 2)
plot!([23,24],[cap_setting(total_residual,0,data["spot"])[end],cap_setting(total_residual,0,data["spot"])[end]],color=:black, label = false, linewidth = 2)

cap_scenario_plot = Plots.plot(highplot, lowplot, constantplot, layout = (1,3), size = (1200,300), bottom_margin = 10*Plots.mm, top_margin=3Plots.mm)
Plots.savefig(cap_scenario_plot,"Figures/cap_scenario.pdf")

### Battery Charging and Discharging Plots ###
charging_discharging_plot = bar(0.5:23.5,transpose(sum(Array(value.(full[0.5,0.5][:pᶜʰ]) .- value.(full[0.5,0.5][:pᵈⁱˢ])) ;dims=1)),bar_width = 1, xlims = (0,24), label = "Battery charging/discharging",ylabel = "Power [kW]",legend =:topleft,color =:black,xticks=0:24)
energy_storage_plot = plot(1:24,transpose(sum(Array(value.(full[0.5,0.5][:e]));dims=1))/sum(data["Eᵐᵃˣ"][1:14]),xlims=(0,24),label = "Battery storage level",legend =:topleft,xlabel = "Time [h]",ylabel = "Battery storage level [%]", color=:black,xticks=0:24)
plot!(0:1,transpose([(transpose(sum(Array(value.(full[0.5,0.5][:e]));dims=1))/sum(data["Eᵐᵃˣ"][1:14]))[end] (transpose(sum(Array(value.(full[0.5,0.5][:e]));dims=1))/sum(data["Eᵐᵃˣ"][1:14]))[1]]),color=:black,label=false)
battery = plot(charging_discharging_plot,energy_storage_plot,layout=(2,1),grid=false)
savefig(battery,"Figures/battery.pdf")

println("       Plotting complete       ")
###################################
### Single run for 14 prosumers ###
###################################

N_consumer = 14 #number of prosumers
grid = "feeder15" #grid to use
data, node, line, ind_cost, tot_cost, ind_profile, tot_profile = data_processing(data_path,grid,N_consumer)

### Residual calculation ###
residual      = sum(eachcol(data["D"] .- data["PV"]))
total_residual  = sum(sum(eachcol(data["D"] .- data["PV"])))

cap_14 = cap_setting(total_residual,1,data["spot"])
community_benchmark, individual_benchmark = benchmarks(ind_cost,ind_profile,tot_profile, cap_14)
b = 0.5

#solving with and without line capacities and checking timing ratio
casestudy14_without_lines = dynamic_pricing(data,node,line,individual_benchmark,b,cap_14,"none",false)
t_14_without = solve_time(casestudy14_without_lines)
obj_14_without = objective_value(casestudy14_without_lines)

casestudy14_with_lines = dynamic_pricing(data,node,line,individual_benchmark,b,cap_14,"none",true)
t_14_with = solve_time(casestudy14_with_lines)
obj_14_with = objective_value(casestudy14_with_lines)

t_ratio = t_14_with/t_14_without

#Line apparent power limit violation check in all hours
k = 0
for t in 1:24
    for l in 1:14
        if value.(casestudy14_without_lines[:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(casestudy14_without_lines[:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 >= (line[l].f̅)^2
            println("Line", i," capacity exceeded in hour", t, "by", value.(casestudy14_without_lines[:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(casestudy14_without_lines[:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 - (line[l].f̅)^2)
            k += 1
        end
    end
end
if k == 0
    println("No Line Violations")
end

###################################
### Single run for 28 prosumers ###
###################################

N_consumer = 28 #number of prosumers
grid = "feeder29" #grid to use
data, node, line, ind_cost, tot_cost, ind_profile, tot_profile = data_processing(data_path,grid,N_consumer)

### Residual calculation ###
residual      = sum(eachcol(data["D"] .- data["PV"]))
total_residual  = sum(sum(eachcol(data["D"] .- data["PV"])))

cap_28 = cap_setting(total_residual,1,data["spot"])
community_benchmark, individual_benchmark = benchmarks(ind_cost,ind_profile,tot_profile, cap_28)

b = 0.5
casestudy28 = dynamic_pricing(data,node,line,individual_benchmark,b,cap_28,"none",false)

t_28 = solve_time(casestudy28)

k = 0
for t in 1:24
    for l in 1:N_consumer
        if value.(casestudy28[:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(casestudy28[:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 >= (line[l].f̅)^2
            println("Line", i," capacity exceeded in hour", t, "by", value.(casestudy28[:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(casestudy28[:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 - (line[l].f̅)^2)
            k += 1
        end
    end
end
if k == 0
    println("No Line Violations")
end

###################################
### Single run for 56 prosumers ###
###################################

N_consumer = 56 #number of prosumers
grid = "feeder57" #grid to use
data, node, line, ind_cost, tot_cost, ind_profile, tot_profile = data_processing(data_path,grid,N_consumer)

### Residual calculation ###
residual      = sum(eachcol(data["D"] .- data["PV"]))
total_residual  = sum(sum(eachcol(data["D"] .- data["PV"])))

cap_56 = cap_setting(total_residual,1,data["spot"])
community_benchmark, individual_benchmark = benchmarks(ind_cost,ind_profile,tot_profile, cap_56)

b = 0.5

casestudy56 = dynamic_pricing(data,node,line,individual_benchmark,b,cap_56,"none",false)
t_56 = solve_time(casestudy56)

for t in 1:24
    for l in 1:N_consumer
        if value.(casestudy56[:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(casestudy56[:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 >= (line[l].f̅)^2
            println("Line", i," capacity exceeded in hour", t, "by", value.(casestudy56[:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(casestudy56[:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 - (line[l].f̅)^2)
            k += 1
        end
    end
end
if k == 0
    println("No Line Violations")
end

###################################
### Single run for 112 prosumers ###
###################################

N_consumer = 112 #number of prosumers
grid = "feeder113" #grid to use
data, node, line, ind_cost, tot_cost, ind_profile, tot_profile = data_processing(data_path,grid,N_consumer)

### Residual calculation ###
residual      = sum(eachcol(data["D"] .- data["PV"]))
total_residual  = sum(sum(eachcol(data["D"] .- data["PV"])))

cap_112 = cap_setting(total_residual,1,data["spot"])
community_benchmark, individual_benchmark = benchmarks(ind_cost,ind_profile,tot_profile, cap_112)

b = 0.5

casestudy112 = dynamic_pricing(data,node,line,individual_benchmark,b,cap_112,"none",false)
t_112 = solve_time(casestudy112)

for t in 1:24
    for l in 1:N_consumer
        if value.(casestudy112[:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(casestudy112[:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 >= (line[l].f̅)^2
            println("Line", i," capacity exceeded in hour", t, "by", value.(casestudy112[:f_p][l,t]/data["Sᵇᵃˢᵉ"]).^2 + value.(casestudy112[:f_q][l,t]/data["Sᵇᵃˢᵉ"]).^2 - (line[l].f̅)^2)
            k += 1
        end
    end
end
if k == 0
    println("No Line Violations")
end

println("14 lines: ",t_14_without)
println("28 lines: ",t_28)
println("56 lines: ",t_56)
println("112 lines: ",t_112)


#Comparative Plots
plot14 = bar(0.5:23.5,value.(casestudy14_without_lines[:pⁱᵐ]),bar_width=1,xlims=(0,24))
bar!(0.5:23.5,value.(casestudy28[:pⁱᵐ]),bar_width=1,xlims=(0,24),color =:green)
plot56 = bar(0.5:23.5,value.(casestudy56[:pⁱᵐ]),bar_width=1,xlims=(0,24))
plot112 = bar(0.5:23.5,value.(casestudy112[:pⁱᵐ]),bar_width=1,xlims=(0,24))

comparison = plot(plot14,plot28,plot56,plot112, layout = (2,2))

bar(0.5:23.5,value.(casestudy112[:pⁱᵐ]),bar_width=1,xlims=(0,24),label="112 Members",xlabel = "Time-of-day [h]",xlabelfontsize=12,xtickfontsize=12,ylabel = "Power [kW]", ylabelfontsize=12,ytickfontsize=12,legendfontsize=12)
bar!(0.5:23.5,value.(casestudy56[:pⁱᵐ]),bar_width=1,xlims=(0,24),color=:turquoise, label="56 Members")
bar!(0.5:23.5,value.(casestudy28[:pⁱᵐ]),bar_width=1,xlims=(0,24),color=:turquoise1, label="28 Members")
bar!(0.5:23.5,value.(casestudy14_without_lines[:pⁱᵐ]),bar_width=1,xlims=(0,24),color=:green, label="14 Members")
plot!(0:23,cap_14,color =:black,linetype=:steppost,linewidth=3,label=false)
plot!(23:24,[cap_14[end],cap_14[end]],color =:black,linetype=:steppost,linewidth=3,label=false)
plot!(0:23,cap_28,color =:black,linetype=:steppost,linewidth=3,label=false)
plot!(23:24,[cap_28[end],cap_14[end]],color =:black,linetype=:steppost,linewidth=3,label=false)
plot!(0:23,cap_56,color =:black,linetype=:steppost,linewidth=3,label=false)
plot!(23:24,[cap_56[end],cap_14[end]],color =:black,linetype=:steppost,linewidth=3,label=false)
plot!(0:23,cap_112,color =:black,linetype=:steppost,linewidth=3,label=false)
plot!(23:24,[cap_112[end],cap_14[end]],color =:black,linetype=:steppost,linewidth=3,label=false)
savefig("Figures/scaling_comparison.pdf")



#Figure generation for latex
plot(1:18, [0.5, 0.5, 0.7, 1, 1.5, 2, 3, 3.5, 3.75, 3.75, 3.5, 3, 2, 1.5, 1, 0.7, 0.5, 0.5],xticks = false, xlabel = "Time", xlabelfontsize = 50, yticks = false, ylabel = "Power", ylabelfontsize=50,grid = false,xlims=(1,18),linewidth=5,linestyle=:dash,label=false,color=:black,left_margin=7Plots.mm)
hline!([3],color =:red, linewidth=5, label = false)
plot!(1:18, [1, 1, 1, 1, 1.5, 2, 3, 3, 3, 3, 3, 3, 2, 1.5, 1, 1, 1, 1], color =:blue, linewidth=5, label = false)

savefig("Figures/illustration.pdf")


plot()
for b in betas
    values_c = []
    for c in caps
        push!(values_c,objective_value(full[b,c]))
    end
    plot!(values_c,label = L"\beta ="*"$b",)
end
display(plot!())

plot()
for c in caps
    values_b = []
    for b in betas
        push!(values_b,objective_value(full[b,c]))
    end
    plot!(values_b,label = L"cap ="*"$c",)
end
display(plot!())