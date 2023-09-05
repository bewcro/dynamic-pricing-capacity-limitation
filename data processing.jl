using CSV, DataFrames, Dates, Query, Plots, DataStructures, ProgressBars

include("DP_CC_OPF/scripts/data_manager.jl")

function dataprocessing_upperlevel(data_path)
        df = CSV.read(data_path*"Load Profile Generator/CREST profiles.csv", DataFrame,header=4,delim=';',skipto=7)
        x = df |> @filter(_.var"Dwelling index" == 1) |> DataFrame

        all_demands = zeros(24,14)
        for j in 1:14
                N = @from i in df begin
                        @where i.var"Dwelling index" == j
                        @select {time = i.Time, demand = i.var"Net dwelling electricity demand"}
                        @collect DataFrame
                end
                N.time = Time.(N.time, "HH.MM.SS p")
                hourly_demand = zeros(24)
                for i in 0:23
                        D = N |> @filter(Hour.(_.time) == Hour(i)) |> DataFrame
                        hourly_demand[i+1] = sum(D[!,"demand"]./(1000*60))
                end
                all_demands[:,j] = hourly_demand
        end
        Demand = all_demands

        solar = CSV.read(string(data_path,"PV.csv"), DataFrame, header = 4)
        solar.local_time = DateTime.(solar.local_time, "yyyy-mm-dd HH:MM")
        pv_base = solar[(solar.local_time .>= DateTime(2019,08,02,00,00,00)) .& (solar."local_time" .< DateTime(2019,08,03,00,00,00)),"electricity"]
        pv = DataFrame(hcat(2*pv_base,zeros(24),zeros(24),pv_base,zeros(24),zeros(24),2*pv_base,2*pv_base,zeros(24),zeros(24),zeros(24),zeros(24),2*pv_base,zeros(24)),:auto)
        prices = DataFrame(CSV.File(string(data_path,"elspotprices.csv")))
        price0208DKK = prices[(prices."HourDK" .>= DateTime(2021,08,02,00,00,00)) .& (prices."HourDK" .< DateTime(2021,08,03,00,00,00)),"SpotPriceDKK"]

        elafgift = 0.7630
        tso = 0.049 + 0.061 + 0.0022
        dso_radius_winter = [0.2296,0.2296,0.2296,0.2296,0.2296,0.2296,0.6889,0.6889,0.6889,0.6889,0.6889,0.6889,0.6889,0.6889,0.6889,0.6889,0.6889,2.0666,2.0666,2.0666,2.0666,0.6889,0.6889,0.6889]
        dso_radius_summer = [0.2296,0.2296,0.2296,0.2296,0.2296,0.2296,0.3444,0.3444,0.3444,0.3444,0.3444,0.3444,0.3444,0.3444,0.3444,0.3444,0.3444,0.8955,0.8955,0.8955,0.8955,0.3444,0.3444,0.3444]
        vindstoed = 0.00375 + 0.000875 + 0.01

        prices = reverse(price0208DKK/1000 .+ elafgift .+ tso) # tariff added https://energinet.dk/El/Elmarkedet/Tariffer/Aktuelle-tariffer
        export_tariff = repeat([vindstoed],24)
        import_tariff = dso_radius_summer

        Delta = 10e-6

        T = 24
        E̅ =    [7.50 10 5.0 5.0 0 10 7.50 5.0 0 10 5.0 5   10 5.0]
        p̅ᶜʰ =  [3.75 5  2.5 2.5 0 5  3.75 2.5 0 5  2.5 2.5 5  2.5]
        p̅ᵈⁱˢ = [3.75 5  2.5 2.5 0 5  3.75 2.5 0 5  2.5 2.5 5  2.5]
        ηᶜʰ = 1.05
        ηᵈⁱˢ = 0.95

        data = Dict("D" => Demand, "PV" => pv, "ηᶜʰ" => ηᵈⁱˢ, "ηᵈⁱˢ" => ηᶜʰ, "αᵍʳⁱᵈ" => 75,#= "βᵍʳⁱᵈ" => 0.5, "Pᵍʳⁱᵈᵐᵃˣ" => p̅ᵍʳⁱᵈ,=#
                    "yⁱᵐ" => import_tariff, "yᵉˣ" => export_tariff, "uᵐᵃˣ" => 1.21, "uᵐⁱⁿ" => 0.81, "spot" => prices, "Eᵐᵃˣ" => E̅, "pᶜʰᵐᵃˣ" => p̅ᶜʰ,
                    "pᵈⁱˢᵐᵃˣ" => p̅ᵈⁱˢ, "M" => 1e4, "pᵇᵃˢᵉ" => 400, "Delta" => Delta)

        caseID = "feeder15"
        (node,line,R,D_n,U_n,U_l,D_l,T) = load_data(caseID)

        ind_cost = Dict()
        tot_cost = 0
        ind_profile = Dict()
        tot_profile = zeros(24)

        residual = sum(eachcol(data["D"] .- data["PV"]))
        total_residual = sum(sum(eachcol(data["D"] .- data["PV"])))
        
        for prosumer_number in 1:length(node)-1

                spot    = data["spot"]
                T       = 24
                PV      = data["PV"][:,prosumer_number]
                Dᵇᵃˢᵉ   = data["D"][:,prosumer_number]
                E̅       = data["Eᵐᵃˣ"][prosumer_number]
                p̅ᶜʰ     = data["pᶜʰᵐᵃˣ"][prosumer_number]
                p̅ᵈⁱˢ    = data["pᵈⁱˢᵐᵃˣ"][prosumer_number]
                ηᶜʰ     = data["ηᶜʰ"]
                ηᵈⁱˢ    = data["ηᵈⁱˢ"]
                yⁱᵐ     = data["yⁱᵐ"]
                yᵉˣ     = data["yᵉˣ"]
                
                opt = optimizer_with_attributes(Gurobi.Optimizer, MOI.Silent() => true)
                prosumer_external = Model(opt)
                @variable(prosumer_external, 0 <= e[1:T] <= E̅)
                @variable(prosumer_external, 0 <= pᶜʰ[1:T] <= p̅ᶜʰ)
                @variable(prosumer_external, 0 <= pᵈⁱˢ[1:T] <= p̅ᵈⁱˢ)
                @variable(prosumer_external, 0 <= p⁺[1:T])
                @variable(prosumer_external, 0 <= p⁻[1:T])
                @variable(prosumer_external, 0 <= q⁺[1:T])
                @variable(prosumer_external, 0 <= q⁻[1:T])
        
                #Objective function
                @objective(prosumer_external, Min, sum((spot[t]+yⁱᵐ[t])*p⁺[t] - (spot[t]-yᵉˣ[t])*p⁻[t] for t in 1:T))
                
                #Constraints
                @constraint(prosumer_external, power_balance[t=1:T], p⁺[t] - p⁻[t] + PV[t] - Dᵇᵃˢᵉ[t] + pᵈⁱˢ[t] - pᶜʰ[t] == 0)
                @constraint(prosumer_external, battery_balance_t1, e[24] + ηᶜʰ*pᶜʰ[1] - ηᵈⁱˢ*pᵈⁱˢ[1] - e[1] == 0)
                @constraint(prosumer_external, battery_balance[t=2:T], e[t-1] + ηᶜʰ*pᶜʰ[t] - ηᵈⁱˢ*pᵈⁱˢ[t] - e[t] == 0)
                @constraint(prosumer_external, active_reactive_plus[t=1:T], q⁺[t] - node[prosumer_number].tan_ϕ*p⁺[t] == 0)
                @constraint(prosumer_external, active_reactive_neg[t=1:T], q⁻[t] - node[prosumer_number].tan_ϕ*p⁻[t] == 0)
                
                set_silent(prosumer_external)
                optimize!(prosumer_external)

                #Cost calculations
                ind_cost[prosumer_number] = objective_value(prosumer_external)
                tot_cost += ind_cost[prosumer_number]
                ind_profile[prosumer_number] = value.(p⁺) - value.(p⁻)
                tot_profile += ind_profile[prosumer_number]
        end
        
        return data, node, line, ind_cost, tot_cost, ind_profile, tot_profile, residual, total_residual
end

function cap_setting(total_residual,var,prices)

        max_price = maximum(prices)
        min_price = minimum(prices)
        mid_price = (max_price + min_price)/2

        avg_cap = total_residual/24

        cap = avg_cap .- var*avg_cap*2*(prices .- mid_price)./(max_price .- min_price)
        return cap
end

function benchmarks(individual_cost,individual_profile,community_profile,cap)
        excess = zeros(length(community_profile))
        for t in 1:length(community_profile)
                excess[t] = max(community_profile[t]-cap[t],0)
        end
        community_penalty = excess .*75
        individual_benchmark = zeros(length(individual_profile))
        for i in 1:length(individual_profile)
                individual_benchmark[i] = individual_cost[i] .+ sum(community_penalty.*(individual_profile[i]./community_profile))
        end

        community_benchmark = sum(individual_benchmark)

        return community_benchmark, individual_benchmark
end

function prosumer_uniform_price(data,node,price,prosumer_number)
    
        spot    = price
        T       = 24
        PV      = data["PV"][:,prosumer_number]
        Dᵇᵃˢᵉ   = data["D"][:,prosumer_number]
        E̅       = data["Eᵐᵃˣ"][prosumer_number]
        p̅ᶜʰ     = data["pᶜʰᵐᵃˣ"][prosumer_number]
        p̅ᵈⁱˢ    = data["pᵈⁱˢᵐᵃˣ"][prosumer_number]
        ηᶜʰ     = data["ηᶜʰ"]
        ηᵈⁱˢ    = data["ηᵈⁱˢ"]
        yⁱᵐ     = data["yⁱᵐ"]
        yᵉˣ     = data["yᵉˣ"]
    
        prosumer_external = Model(Gurobi.Optimizer)
        @variable(prosumer_external, 0 <= e[t=1:T] <= E̅)
        @variable(prosumer_external, 0 <= pᶜʰ[t=1:T] <= p̅ᶜʰ)
        @variable(prosumer_external, 0 <= pᵈⁱˢ[t=1:T] <= p̅ᵈⁱˢ)

        @variable(prosumer_external, 0 <= p⁺[1:T])
        @variable(prosumer_external, 0 <= p⁻[1:T])
        @variable(prosumer_external, 0 <= q⁺[1:T])
        @variable(prosumer_external, 0 <= q⁻[1:T])
    
        #Objective function
        @objective(prosumer_external, Min, sum((spot[t]+yⁱᵐ[t])*p⁺[t] - (spot[t]-yᵉˣ[t])*p⁻[t] for t in 1:T))
        
        #Constraints
        @constraint(prosumer_external, power_balance[t=1:T], p⁺[t] - p⁻[t] + PV[t] - Dᵇᵃˢᵉ[t] + + pᵈⁱˢ[t] - pᶜʰ[t] == 0)
        @constraint(prosumer_external, battery_balance_t1[t=1], e[24] + ηᶜʰ*pᶜʰ[t] - ηᵈⁱˢ*pᵈⁱˢ[t] - e[t] == 0)
        @constraint(prosumer_external, battery_balance[t=2:T], e[t-1] + ηᶜʰ*pᶜʰ[t] - ηᵈⁱˢ*pᵈⁱˢ[t] - e[t]== 0)
        @constraint(prosumer_external, active_reactive_plus[t=1:T], q⁺[t] - node[prosumer_number].tan_ϕ*p⁺[t] == 0)
        @constraint(prosumer_external, active_reactive_neg[t=1:T], q⁻[t] - node[prosumer_number].tan_ϕ*p⁻[t] == 0)
        
        optimize!(prosumer_external)
        return prosumer_external
    end
