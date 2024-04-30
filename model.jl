using JuMP, Gurobi

function dynamic_pricing(data, node, line, individual, beta, cap, distribution, line_cap)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 1e-4)
    set_optimizer_attribute(model, "OutputFlag", 0)
    # Sets
    I = 1:length(node)-1
    N = length(node)-1
    T = 24
    ##################
    ### PARAMETERS ###
    ##################

    Dᵇᵃˢᵉ   = data["D"]         #Demand [kWh]
    PV      = data["PV"]        #Solar power [kWh]
    η       = data["η"]         #Charge/Discharge efficiency [%]
    αᵍʳⁱᵈ   = data["αᵍʳⁱᵈ"]     #penalty for exceeding capacity limit [DKK]
    yⁱᵐ     = data["yⁱᵐ"]       #DSO import tariff [DKK]
    yᵉˣ     = data["yᵉˣ"]       #DSO export tariff [DKK]
    u̅       = data["uᵐᵃˣ"]      #maximum allowed voltage [p.u]
    u̲       = data["uᵐⁱⁿ"]      #minimum allowed voltage [p.u]
    spot    = data["spot"]      #DA spot prices [DKK]
    E̅       = data["Eᵐᵃˣ"]      #Battery capacity [kWh]
    p̅ᶜʰ     = data["pᶜʰᵐᵃˣ"]    #Maximum battery charge power [kW]
    p̅ᵈⁱˢ    = data["pᵈⁱˢᵐᵃˣ"]   #Maximum battery discharge power [kW]
    Sᵇᵃˢᵉ   = data["Sᵇᵃˢᵉ"]
    γ       = data["Delta"]     #Penalization factor regularisation factor [no unit]
    βᵍʳⁱᵈ = beta                #Discount factor for internal power flows [%]
    VOLL =  αᵍʳⁱᵈ*1.25               #Value of lost load (> grid penalty) [DKK]
    P̅ᵍʳⁱᵈ = cap                 #Capacity limit imposed by DSO [kW]

    #################
    ### VARIABLES ###
    #################
    @variable(model, 0 <= x[I,1:T])                     #price [DKK]
    @variable(model, x̅)                                 #maximum price [DKK]
    @variable(model, 0 <= pⁱᵐ[1:T])                     #power imported from grid [kW]
    @variable(model, 0 <= pᵉˣ[1:T])                     #power exported to grid [kW]
    @variable(model, 0 <= qⁱᵐ[1:T])                     #reactive power imported by community [kW]
    @variable(model, 0 <= qᵉˣ[1:T])                     #reactive power exported [kW]
    @variable(model, 0 <= pᵖᵉⁿ[1:T])                    #amount of power above capacity limit [kW]
    @variable(model, 0 <= e[i=I,t=1:T] <= E̅[i])         #battery charge level [kWh]
    @variable(model, 0 <= pᶜʰ[i=I,t=1:T] <= p̅ᶜʰ[i])     #power charged to battery [kW]
    @variable(model, 0 <= pᵈⁱˢ[i=I,t=1:T] <= p̅ᵈⁱˢ[i])   #power discharged from battery [kW]
    @variable(model, 0 <= p⁺[I,1:T])                    #power imported by prosumer i [kW]
    @variable(model, 0 <= p⁻[I,1:T])                    #power exported by prosumer i [kW]
    @variable(model, 0 <= q⁺[I,1:T])                    #reactive power import per prosumer [kVAr]
    @variable(model, 0 <= q⁻[I,1:T])                    #reactice power export per prosumer
    @variable(model, f_p[1:N,1:T])                      #active power flow downstream from node n [kW]
    @variable(model, f_q[1:N,1:T])                      #reactive power flow downstream from node n [kW]
    @variable(model, u[0:N,1:T])                        #square voltage at node n [p.u]
    @variable(model, 0 <= dˢʰᵉᵈ[I,1:T])                 #load shedding at node n [kW]

    #equality constraint duals
    @variable(model, λ₁[I,1:T])                         #dual for power balance
    @variable(model, λ₂[I])                             #dual for battery balance time step 1
    @variable(model, λ₃[I,2:T])                         #dual for battery balance
    @variable(model, λ₄[I,1:T])                         #dual for active/reactive power import equality
    @variable(model, λ₅[I,1:T])                         #dual for active/reactive power export equality

    #inequality constraint duals
    @variable(model, μ₁[I,1:T] >= 0)                    #lower bound p⁺
    @variable(model, μ₂[I,1:T] >= 0)                    #lower bound p⁻
    @variable(model, μ₃[I,1:T] >= 0)                    #lower bound q⁺
    @variable(model, μ₄[I,1:T] >= 0)                    #lower bound q⁻
    @variable(model, μ₅[I,1:T] >= 0)                    #lower bound pᶜʰ
    @variable(model, μ₆[I,1:T] >= 0)                    #lower bound pᵈⁱˢ
    @variable(model, μ₇[I,1:T] >= 0)                    #lower bound E
    @variable(model, μ₈[I,1:T] >= 0)                    #lower bound dˢʰᵉᵈ
    @variable(model, μ₉[I,1:T] >= 0)                    #upper bound pᶜʰ
    @variable(model, μ₁₀[I,1:T] >= 0)                   #upper bound pᵈⁱˢ
    @variable(model, μ₁₁[I,1:T] >= 0)                   #upper bound E
    @variable(model, μ₁₂[I,1:T] >= 0)                   #upper bound dˢʰᵉᵈ

    #SOS1 complementarity linearization variables
    @variable(model, u₁[I,1:T])                         #lower bound p⁺
    @variable(model, u₂[I,1:T])                         #lower bound p⁻
    @variable(model, u₃[I,1:T])                         #lower bound q⁺
    @variable(model, u₄[I,1:T])                         #lower bound q⁻
    @variable(model, u₅[I,1:T])                         #lower bound pᶜʰ
    @variable(model, u₆[I,1:T])                         #lower bound pᵈⁱˢ
    @variable(model, u₇[I,1:T])                         #lower bound E
    @variable(model, u₈[I,1:T])                         #lower bound dˢʰᵉᵈ
    @variable(model, u₉[I,1:T])                         #upper bound pᶜʰ
    @variable(model, u₁₀[I,1:T])                        #upper bound pᵈⁱˢ
    @variable(model, u₁₁[I,1:T])                        #upper bound E
    @variable(model, u₁₂[I,1:T])                        #upper bound dˢʰᵉᵈ
    @variable(model, v₁[I,1:T,1:2])
    @constraint(model, v₁_sos1[i=I,t=1:T], [v₁[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₂[I,1:T,1:2])
    @constraint(model, v₂_sos1[i=I,t=1:T], [v₂[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₃[I,1:T,1:2])
    @constraint(model, v₃_sos1[i=I,t=1:T], [v₃[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₄[I,1:T,1:2])
    @constraint(model, v₄_sos1[i=I,t=1:T], [v₄[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₅[I,1:T,1:2])
    @constraint(model, v₅_sos1[i=I,t=1:T], [v₅[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₆[I,1:T,1:2])
    @constraint(model, v₆_sos1[i=I,t=1:T], [v₆[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₇[I,1:T,1:2])
    @constraint(model, v₇_sos1[i=I,t=1:T], [v₇[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₈[I,1:T,1:2])
    @constraint(model, v₈_sos1[i=I,t=1:T], [v₈[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₉[I,1:T,1:2])
    @constraint(model, v₉_sos1[i=I,t=1:T], [v₉[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₁₀[I,1:T,1:2])
    @constraint(model, v₁₀_sos1[i=I,t=1:T], [v₁₀[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₁₁[I,1:T,1:2])
    @constraint(model, v₁₁_sos1[i=I,t=1:T], [v₁₁[i,t,s] for s in 1:2] in SOS1())
    @variable(model, v₁₂[I,1:T,1:2])
    @constraint(model, v₁₂_sos1[i=I,t=1:T], [v₁₁[i,t,s] for s in 1:2] in SOS1())

    #individual rationality variables
    @variable(model, ω⁺[I] >= 0)                        #slack variable for positive deviation [DKK]
    @variable(model, ω⁻[I] >= 0)                        #slack variable for negative deviation [DKK]
    @variable(model, Ω[1:2] >= 0)                       #Variable to ensure only positive or negative slack

    #Constraints
    #individual rationality constraints
    @constraint(model, individual_rationality[i=I], sum(λ₁[i,t]*(PV[t,i]-Dᵇᵃˢᵉ[t,i]) - μ₉[i,t]*p̅ᶜʰ[i] - μ₁₀[i,t]*p̅ᵈⁱˢ[i] - μ₁₁[i,t]*E̅[i] - μ₁₂[i,t]*Dᵇᵃˢᵉ[t,i] - VOLL*dˢʰᵉᵈ[i,t] for t=1:T) == individual[i] + ω⁺[i] - ω⁻[i])
    @constraint(model, slack_sos1, [Ω[s] for s in 1:2] in SOS1())
    @constraint(model, slack_equality_plus, sum(ω⁺[i] for i=I) == Ω[1])
    @constraint(model, slack_equality_minus, sum(ω⁻[i] for i=I) == Ω[2])

    #Power penalty constraint
    @constraint(model, power_penalty[t=1:T], pᵖᵉⁿ[t] >= pⁱᵐ[t] - P̅ᵍʳⁱᵈ[t])

    #Cost balance constraint
    @constraint(model, cost_balance, sum(λ₁[i,t]*(PV[t,i]-Dᵇᵃˢᵉ[t,i]) - μ₉[i,t]*p̅ᶜʰ[i] - μ₁₀[i,t]*p̅ᵈⁱˢ[i] - μ₁₁[i,t]*E̅[i] - μ₁₂[i,t]*Dᵇᵃˢᵉ[t,i] - VOLL*dˢʰᵉᵈ[i,t] for t=1:T for i in I) == sum(pⁱᵐ[t]*(spot[t] + yⁱᵐ[t]) - pᵉˣ[t]*(spot[t] - yᵉˣ[t]) + αᵍʳⁱᵈ*pᵖᵉⁿ[t] + (1-βᵍʳⁱᵈ)*yⁱᵐ[t]*(sum(p⁺[i,t] for i=I) - pⁱᵐ[t]) for t=1:T))

    #price regularisation constraint
    @constraint(model, price_regularisation[i=I,t=1:T], x[i,t] <= x̅)
    
    #LinDistFlow constraints - cite M.Baran - 
    @constraint(model, activecommunitybalance[n=0,t=1:T], pⁱᵐ[t] - pᵉˣ[t] == sum(f_p[j,t] for j in node[n].C))
    @constraint(model, reactivecommunitybalance[n=0,t=1:T], qⁱᵐ[t] - qᵉˣ[t] == sum(f_q[j,t] for j in node[n].C))
    @constraint(model, activepower[i=I, t=1:T], f_p[i,t] == p⁺[i,t] - p⁻[i,t] + sum(f_p[j,t] for j in node[i].C))
    @constraint(model, reactivepower[i=I, t=1:T], f_q[i,t] == q⁺[i,t] - q⁻[i,t] + sum(f_q[j,t] for j in node[i].C))
    @constraint(model, voltage[n=1:N, t=1:T], u[node[n].A[1],t] - 2*((f_p[n,t]/Sᵇᵃˢᵉ)*line[n].r + (f_q[n,t]/Sᵇᵃˢᵉ)*line[n].x) == u[n,t])
    @constraint(model, base_voltage[t=1:T], u[0,t] == 1)
    @constraint(model, voltagemax[n=0:N, t=1:T], u̲ <= u[n,t] <= u̅)
    #line flow constraint is eliminated to maintain LinDistFlow formulation and ease computational burden
    if line_cap == true
        @constraint(model, linelimit[l=1:N, t=1:T], (f_p[l,t]/Sᵇᵃˢᵉ)^2 + (f_q[l,t]/Sᵇᵃˢᵉ)^2 <= line[l].f̅^2)
    end

    #Karush Kuhn Tucker Optimality Conditions
    @constraint(model, KKT_pplus[i=I,t=1:T], x[i,t] + λ₁[i,t] - node[i].tan_ϕ*λ₄[i,t] - μ₁[i,t] == 0)
    @constraint(model, KKT_pminus[i=I,t=1:T], -x[i,t] - λ₁[i,t] - node[i].tan_ϕ*λ₅[i,t] - μ₂[i,t] == 0)
    @constraint(model, KKT_qplus[i=I,t=1:T], λ₄[i,t] - μ₃[i,t] == 0)
    @constraint(model, KKT_qminus[i=I,t=1:T], λ₅[i,t] - μ₄[i,t] == 0)
    @constraint(model, KKT_pch1[i=I], -λ₁[i,1]/η + λ₂[i] - μ₅[i,1] + μ₉[i,1] == 0)
    @constraint(model, KKT_pch[i=I,t=2:T], -λ₁[i,t]/η + λ₃[i,t] - μ₅[i,t] + μ₉[i,t] == 0)
    @constraint(model, KKT_pdis1[i=I], η*λ₁[i,1] - λ₂[i] - μ₆[i,1] + μ₁₀[i,1] == 0)
    @constraint(model, KKT_pdis[i=I,t=2:T], η*λ₁[i,t] - λ₃[i,t] - μ₆[i,t] + μ₁₀[i,t] == 0)
    @constraint(model, KKT_E1[i=I], -λ₂[i] + λ₃[i,2] - μ₇[i,1] + μ₁₁[i,1] == 0)
    @constraint(model, KKT_Et[i=I,t=2:T-1], -λ₃[i,t] + λ₃[i,t+1] - μ₇[i,t] + μ₁₁[i,t] == 0)
    @constraint(model, KKT_Etmax[i=I], -λ₃[i,T] + λ₂[i] - μ₇[i,T] + μ₁₁[i,T] == 0)
    @constraint(model, KKT_dshed[i=I,t=1:T], VOLL + λ₁[i,t] - μ₈[i,t] + μ₁₂[i,t] == 0)
    @constraint(model, KKT_bat_λ₁[i=I,t=1:T], p⁺[i,t] - p⁻[i,t] + PV[t,i] - Dᵇᵃˢᵉ[t,i] + dˢʰᵉᵈ[i,t] + η*pᵈⁱˢ[i,t] - pᶜʰ[i,t]/η == 0)
    @constraint(model, KKT_bat_λ₂[i=I],   e[i,24] + pᶜʰ[i,1] - pᵈⁱˢ[i,1] - e[i,1] == 0)
    @constraint(model, KKT_bat_λ₃[i=I,t=2:T], e[i,t-1] + pᶜʰ[i,t] - pᵈⁱˢ[i,t] - e[i,t]== 0)
    @constraint(model, KKT_bat_λ₄[i=I,t=1:T], q⁺[i,t] - node[i].tan_ϕ*p⁺[i,t] == 0)
    @constraint(model, KKT_bat_λ₅[i=I,t=1:T], q⁻[i,t] - node[i].tan_ϕ*p⁻[i,t] == 0)

    #SOS1 reformulation
    @constraint(model, complementarity_μ₁_1[i=I,t=1:T], u₁[i,t] == (μ₁[i,t] + p⁺[i,t])/2)
    @constraint(model, complementarity_μ₁_2[i=I,t=1:T], v₁[i,t,1] - v₁[i,t,2] == (μ₁[i,t] - p⁺[i,t])/2)
    @constraint(model, complementarity_μ₁_3[i=I,t=1:T], u₁[i,t] - (v₁[i,t,1] + v₁[i,t,2]) == 0)

    @constraint(model, complementarity_μ₂_1[i=I,t=1:T], u₂[i,t] == (μ₂[i,t] + p⁻[i,t])/2)
    @constraint(model, complementarity_μ₂_2[i=I,t=1:T], v₂[i,t,1] - v₂[i,t,2] == (μ₂[i,t] - p⁻[i,t] )/2)
    @constraint(model, complementarity_μ₂_3[i=I,t=1:T], u₂[i,t] - (v₂[i,t,1] + v₂[i,t,2]) == 0)

    @constraint(model, complementarity_μ₃_1[i=I,t=1:T], u₃[i,t] == (μ₃[i,t] + q⁺[i,t])/2)
    @constraint(model, complementarity_μ₃_2[i=I,t=1:T], v₃[i,t,1] - v₃[i,t,2] == (μ₃[i,t] - q⁺[i,t])/2)
    @constraint(model, complementarity_μ₃_3[i=I,t=1:T], u₃[i,t] - (v₃[i,t,1] + v₃[i,t,2]) == 0)

    @constraint(model, complementarity_μ₄_1[i=I,t=1:T], u₄[i,t] == (μ₄[i,t] + q⁻[i,t])/2)
    @constraint(model, complementarity_μ₄_2[i=I,t=1:T], v₄[i,t,1] - v₄[i,t,2] == (μ₄[i,t] - q⁻[i,t])/2)
    @constraint(model, complementarity_μ₄_3[i=I,t=1:T], u₄[i,t] - (v₄[i,t,1] + v₄[i,t,2]) == 0)
    
    @constraint(model, complementarity_μ₅_1[i=I,t=1:T], u₅[i,t] == (μ₅[i,t] + pᶜʰ[i,t])/2)
    @constraint(model, complementarity_μ₅_2[i=I,t=1:T], v₅[i,t,1] - v₅[i,t,2] == (μ₅[i,t] - pᶜʰ[i,t])/2)
    @constraint(model, complementarity_μ₅_3[i=I,t=1:T], u₅[i,t] - (v₅[i,t,1] + v₅[i,t,2]) == 0)

    @constraint(model, complementarity_μ₆_1[i=I,t=1:T], u₆[i,t] == (μ₆[i,t] + pᵈⁱˢ[i,t])/2)
    @constraint(model, complementarity_μ₆_2[i=I,t=1:T], v₆[i,t,1] - v₆[i,t,2] == (μ₆[i,t] - pᵈⁱˢ[i,t])/2)
    @constraint(model, complementarity_μ₆_3[i=I,t=1:T], u₆[i,t] - (v₆[i,t,1] + v₆[i,t,2]) == 0)

    @constraint(model, complementarity_μ₇_1[i=I,t=1:T], u₇[i,t] == (μ₇[i,t] + e[i,t])/2)
    @constraint(model, complementarity_μ₇_2[i=I,t=1:T], v₇[i,t,1] - v₇[i,t,2] == (μ₇[i,t] - e[i,t])/2)
    @constraint(model, complementarity_μ₇_3[i=I,t=1:T], u₇[i,t] - (v₇[i,t,1] + v₇[i,t,2]) == 0)

    @constraint(model, complementarity_μ₈_1[i=I,t=1:T], u₈[i,t] == (μ₈[i,t] + dˢʰᵉᵈ[i,t])/2)
    @constraint(model, complementarity_μ₈_2[i=I,t=1:T], v₈[i,t,1] - v₈[i,t,2] == (μ₈[i,t] - dˢʰᵉᵈ[i,t])/2)
    @constraint(model, complementarity_μ₈_3[i=I,t=1:T], u₈[i,t] - (v₈[i,t,1] + v₈[i,t,2]) == 0)

    @constraint(model, complementarity_μ₉_1[i=I,t=1:T], u₉[i,t] == (μ₉[i,t] + (p̅ᶜʰ[i] - pᶜʰ[i,t]))/2)
    @constraint(model, complementarity_μ₉_2[i=I,t=1:T], v₉[i,t,1] - v₉[i,t,2] == (μ₉[i,t] - (p̅ᶜʰ[i] - pᶜʰ[i,t]))/2)
    @constraint(model, complementarity_μ₉_3[i=I,t=1:T], u₉[i,t] - (v₉[i,t,1] + v₉[i,t,2]) == 0)

    @constraint(model, complementarity_μ₁₀_1[i=I,t=1:T], u₁₀[i,t] == (μ₁₀[i,t] + (p̅ᵈⁱˢ[i] - pᵈⁱˢ[i,t]))/2)
    @constraint(model, complementarity_μ₁₀_2[i=I,t=1:T], v₁₀[i,t,1] - v₁₀[i,t,2] == (μ₁₀[i,t] - (p̅ᵈⁱˢ[i] - pᵈⁱˢ[i,t]))/2)
    @constraint(model, complementarity_μ₁₀_3[i=I,t=1:T], u₁₀[i,t] - (v₁₀[i,t,1] + v₁₀[i,t,2]) == 0)

    @constraint(model, complementarity_μ₁₁_1[i=I,t=1:T], u₁₁[i,t] == (μ₁₁[i,t] + (E̅[i] - e[i,t]))/2)
    @constraint(model, complementarity_μ₁₁_2[i=I,t=1:T], v₁₁[i,t,1] - v₁₁[i,t,2] == (μ₁₁[i,t] - (E̅[i] - e[i,t]))/2)
    @constraint(model, complementarity_μ₁₁_3[i=I,t=1:T], u₁₁[i,t] - (v₁₁[i,t,1] + v₁₁[i,t,2]) == 0)

    @constraint(model, complementarity_μ₁₂_1[i=I,t=1:T], u₁₂[i,t] == (μ₁₂[i,t] + (Dᵇᵃˢᵉ[t,i] - dˢʰᵉᵈ[i,t]))/2)
    @constraint(model, complementarity_μ₁₂_2[i=I,t=1:T], v₁₂[i,t,1] - v₁₂[i,t,2] == (μ₁₂[i,t] - (Dᵇᵃˢᵉ[t,i] - dˢʰᵉᵈ[i,t]))/2)
    @constraint(model, complementarity_μ₁₂_3[i=I,t=1:T], u₁₂[i,t] - (v₁₂[i,t,1] + v₁₂[i,t,2]) == 0)

    #Objective function definition depending on which is selected
    if distribution == "none"
        #Objective Function
        @objective(model, Min, sum(pⁱᵐ[t]*(spot[t] + yⁱᵐ[t]) - pᵉˣ[t]*(spot[t] - yᵉˣ[t]) + αᵍʳⁱᵈ*pᵖᵉⁿ[t] + (1-βᵍʳⁱᵈ)*yⁱᵐ[t]*(sum(p⁺[i,t] for i=I) - pⁱᵐ[t]) for t=1:T) + γ*x̅^2 + sum(VOLL*dˢʰᵉᵈ[i,t] for i=I for t=1:T) + 100*sum(ω⁺[i] for i=I))
    elseif distribution == "uniform"
        #Objective Function
        @objective(model, Min, sum(pⁱᵐ[t]*(spot[t] + yⁱᵐ[t]) - pᵉˣ[t]*(spot[t] - yᵉˣ[t]) + αᵍʳⁱᵈ*pᵖᵉⁿ[t] + (1-βᵍʳⁱᵈ)*yⁱᵐ[t]*(sum(p⁺[i,t] for i=I) - pⁱᵐ[t]) for t=1:T) + γ*x̅^2 + sum(VOLL*dˢʰᵉᵈ[i,t] for i=I for t=1:T) + 100*sum(ω⁺[i] for i=I) + γ*sum((ω⁺[i] - sum(ω⁺[i] for i=I)/length(I) + ω⁻[i] - sum(ω⁻[i] for i=I)/length(I))^2 for i=I))
    elseif distribution == "proportional"
        #Objective Function
        @objective(model, Min, sum(pⁱᵐ[t]*(spot[t] + yⁱᵐ[t]) - pᵉˣ[t]*(spot[t] - yᵉˣ[t]) + αᵍʳⁱᵈ*pᵖᵉⁿ[t] + (1-βᵍʳⁱᵈ)*yⁱᵐ[t]*(sum(p⁺[i,t] for i=I) - pⁱᵐ[t]) for t=1:T) + γ*x̅^2 + sum(VOLL*dˢʰᵉᵈ[i,t] for i=I for t=1:T) + 100*sum(ω⁺[i] for i=I) + γ*sum((ω⁺[i] + ω⁻[i] - (sum(ω⁺[i] + ω⁻[i] for i=I)*sum(Dᵇᵃˢᵉ[t,i] - PV[t,i] for t in 1:T))/(sum(Dᵇᵃˢᵉ[t,i] - PV[t,i] for i in I, t in 1:T)))^2 for i=I))
    end

    optimize!(model)
    println("--------------------")
    println(objective_value(model))
    println("--------------------")
    return model
end 