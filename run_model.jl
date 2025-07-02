#!/usr/bin/env julia

using Plots
using Printf
using DataStructures: OrderedDict
using NCDatasets
using Arrow, DataFrames
using CSV, DataFrames
using Colors
using ColorSchemes
using Plots
using Printf
using LaTeXStrings
using Profile
using Statistics 

include("calculate_physical_variables.jl") 
include("advance_variables.jl")
include("phytoplankton.jl")
include("forcings.jl") 
include("output.jl")
include("calculate_heat_fluxes.jl")

# ws1 = parse(Float64, ARGS[1])
# ws2 = parse(Float64, ARGS[2])
# pmax1 = parse(Float64, ARGS[3])
# pmax2 = parse(Float64, ARGS[4])
# fout_name = ARGS[5]

function run_my_model(ws1::Real, ws2::Real, pmax1::Real, pmax2::Real, file_out_name::String)

    println("Running model with ws1 = $ws1, ws2 = $ws2, pmax1 = $pmax1, pmax2 = $pmax2.. \n output file name = $file_out_name \n")

    forcing_fn = "/global/homes/s/siennaw/scratch/siennaw/stockton_field_data/forcing_for_model/august6-28/"
    # forcing_fn = "/global/homes/s/siennaw/scratch/siennaw/stockton_field_data/forcing_for_model/2022/june2-29"
    #***********************************************************************
    # Read in the feather file 
    # data  = Arrow.Table("interpolated_temperature_profile_aug10-16.feather")
    # temp_data = DataFrame(data)

    # function get_temp_field(index::Int, temp_data=temp_data)
    #     return collect(temp_data[xf])
    # end

    # function get_unstrat_temp_field(index::Int, temp_data=temp_data)
    #     return fill( mean(collect(temp_data[index,1:end-1])), 60) 
    # end

    #***********************************************************************
    # Read in the heat flux data 
    # heat_flux_fn =  @sprintf("%s/heat_flux.csv", forcing_fn)  

    forcing_fp = @sprintf("%s/CIMIS.csv", forcing_fn)  
    @info ("Reading in forcing data from $forcing_fp ...")
    forcing_df = CSV.read(forcing_fp, DataFrame)
    

    heat_flux_fn =  @sprintf("%s/heat_flux2.csv", forcing_fn)  
    @info ("Reading in heat flux data from $heat_flux_fn ...")
    dfhf = CSV.read(heat_flux_fn, DataFrame)
    # surface0 = df[!,"surface"]
    # shortwave = df[!,"sw_in"]

    function get_surf_flux(index::Int, dfhf=dfhf)
        sf = dfhf[index,"surface"]/(specific_heat_water * rhoW) 
        return sf
    end

    function get_shortwave(index::Int, dfhf=dfhf)
        sw = dfhf[index,"sw_in"]
        return sw 
    end

    #***********************************************************************
    # Read in the CIMIS data
    cimis_fn =  @sprintf("%s/PAR.csv", forcing_fn)  
    @info ("Reading in CIMIS data from $heat_flux_fn ...")

    df = CSV.read(cimis_fn, DataFrame)
    par = df[!,"Sol Rad (PAR)"]

    function get_light(index::Int, par=par)
        return par[index]
    end
    #***********************************************************************


    #***********************************************************************
    # Wind time series 
    wind_fn =  @sprintf("%s/wind.csv", forcing_fn)  
    @info ("Reading in wind data from $wind_fn ...")
    df = CSV.read(wind_fn, DataFrame)
    wind = df[!,"WindSpeed"] 
    real_time = df[!,"time"]
    
    
    function get_wind_speed(index::Int, wind=wind)
        return wind[index]
    end
    #***********************************************************************


    #********************** SPATIAL DOMAIN  ***************************
    N = 60    # number of grid points
    H = 6    # depth (meters)
    dz = H/N  # grid spacing - may need to adjust to reduce oscillations
    dt = 10   # (seconds) size of time step 
    M  = 3600*16 #359*24*22 #206000 # 359*24* 25 # 17280 #00 #000 # 50000  #500 # 17280 # 
         
    # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
    isave = 100 # 360 #1000
    var2save = ["U","Kq", "Nu", "C", "Kz", "L", "Q2", "Q2L", "N_BV2", "algae1","algae2"]

    create_output_dict(M, isave, var2save, N)

    # Create depth vector 
    z = collect(H:-dz:dz) .- dz/2 # depth vector

    function get_body_flux_heat(sw, N, z)
        rad = zeros(N)       # radiation 
        qsource = zeros(N)   # heat source term
        # A = 0.78             # From Paulson & Simpson for Type III water 
        # g1 = 1.4
        # g2 = 7.9

        # try 
        g1 = 0.33
        g2 = 2.34
        A = 0.57
        dz = 0.1
        rhoW = 1000                     # Density of water, kg/m^3
        specific_heat_water = 4181 
        # sw = sw * 1.3
        for i in 1:N 
            rad[i] = sw * (A * exp(-z[i]/g1) + (1-A)*exp(-z[i]/g2))
        end 
        qsource[end] = (sw - rad[end - 1])/(specific_heat_water * rhoW * dz)  # dz  z[end] 
        for i in range(2, N)
            qsource[i] = (rad[i] - rad[i - 1])/(specific_heat_water * rhoW * dz) # mark thinks this should be dz 
        end 
        return qsource
        
    end

    #********************** FIXED CONSTANTS  ***************************
    rhoA = 1.23                     # Density of air, kg/m^3
    rhoW = 1000                     # Density of water, kg/m^3
    specific_heat_water = 4181      # J/kg-degC
    specific_heat_air = 1007        # J/kg-degC x RH
    c_d = 0.05                      # Drag coefficient [-]
    cm2m = 0.01
    hr2s = 1/3600

    #********************** INITIAL CONDITION ***************************
    # Initialize thermocline based on tanh curve 
    base_temp = 28
    dtemp = 1.5 
    stretch = 0.25 

    #********************** DEFINE HYDRODYNAMIC FORCINGS ***************************
    # (1) PRESSURE 
    Px0 = 2e-6          # Pressure gradient forcing
    T_Px = 12           # Period [hours] on pressure gradient forcing. Set to 0 for steady

    # (2) Wind
    # Wind = 1                       # u_star =m/s >> 0.05 is  drag coefficient, 10 is my wind speed 
    # WIND = (c_d * Wind)^2 * rhoA   # this is rho * u*^2

    # (3) Temperature
    top_temp = 28.2
    bottom_temp = 27.8

    # (4) Light 
    DIURNAL_LIGHT = false  
    background_turbidity =  0.5
    # I_in = 350 


    #********************** DEFINE PHYTOPLANKTON FORCINGS ***************************
    # init_algae = 0.005
    init_algae = 0.01

    algae1 = Dict("k" => 0.034,              # specific light attenuation coefficient [cm^2 / 10^6 cells]
    "pmax" =>pmax1 * hr2s,           # maximum specific growth rate [1/hour]
    "ws" => 1.38e-4, #1.38e-4,           # vertical velocity [m/s]
    "Hi" => 40,                # half-saturation of light-limited growth [mu mol photons * m^2/s]
    "Li" => 0.004 * hr2s,             # specific loss rate [1/hour]
    "name" => "HAB",           # name of the species
    "self_shading" => true)    # self-shading effect (true/false)


    # algae1 = Dict("k" => 0.7,              # specific light attenuation coefficient [cm^2 / 10^6 cells]
    #             "pmax" => pmax1, #0.05 * hr2s,          # maximum specific growth rate [1/hour]
    #             "ws" => ws1, #-1.38e-5,           # vertical velocity [m/s]
    #             "Hi" => 40,                 # half-saturation of light-limited growth [mu mol photons * m^2/s]
    #             "Li" => 0.006 * hr2s,             # specific loss rate [1/hour]
    #             "name" => "Diatom",           # name of the species
    #             "self_shading" => true)    # self-shading effect (true/false))       

    algae2 = Dict("k" => 0.034,              # specific light attenuation coefficient [cm^2 / 10^6 cells]
                "pmax" => pmax2, #0.008 * hr2s,           # maximum specific growth rate [1/hour]
                "ws" => ws2, #1.38e-4,           # vertical velocity [m/s]
                "Hi" => 40,                # half-saturation of light-limited growth [mu mol photons * m^2/s]
                "Li" => 0.004 * hr2s,             # specific loss rate [1/hour]
                "name" => "HAB",           # name of the species
                "self_shading" => true)    # self-shading effect (true/false)
    # '''
    # Diatoms ws = -1.38e-5 m/s
    # Cyanobacteria = 1.38e-4 m/s
    # '''
    # get_unstrat_temp_field(1) # 
    #***************************************************************************
    #   Initialize variables
    #***************************************************************************

    # Create dictionary to hold important discretization parameters
    discretization = Dict("beta" => (dt/dz^2), "dz" => dz, "dt" => dt, "N" => N, "z"=> z, "H" => H)

    # Initalize velocity
    U = similar(z) .+ 1e-1
    # C = zeros(N) .+ LinRange(bottom_temp, top_temp, N)  
    # C = @. base_temp - tanh(z * 2)*0.4 

    C = [ 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.798220, 27.816607, 27.834993, 27.853380, 27.878350, 27.915170, 27.951990, 27.988810, 28.025630, 28.062450, 28.099270, 28.136090, 28.172910, 28.181460, 28.182300, 28.183140, 28.183980, 28.184820, 28.185660, 28.186500, 28.187340, 28.188050, 28.188283, 28.188517, 28.188750, 28.188983, 28.189087, 28.189180, 28.189273, 28.189367, 28.189280, 28.189093, 28.188907, 28.188720, 28.188400, 28.187840, 28.187280, 28.186720, 28.186233, 28.186700, 28.187167, 28.187027, 28.186840, 28.188853, 28.191467, 28.189520, 28.185040,]
    # C = @. C/max(C) * 26
    # C = zeros(N) .+ 26
    rho = calculate_rho(C, base_temp) 
    N_BV2 = calculate_brunt_vaisala(rho, discretization)

    algae1["c"] = zeros(N) .+ init_algae 
    algae2["c"] = zeros(N) .+ init_algae 

    Q2, Q2L, Q, L, Gh, nu_t, Kq, Kz = initialize_turbulent_functions(discretization, N_BV2)

    # Initial dictionary to store variables
    variables = Dict("U" => U, "C" => C, "N_BV2" => N_BV2, 
                    "Nu" => nu_t, "Q2" => Q2, "Q2L" => Q2L, 
                    "Kq" => Kq, "Kz" => Kz, "L" => L)

    Times = collect(0:dt:(M*dt))
    real_times_saved = []

    #***************************************************************************
    save2output(1, 1, "algae1", algae1["c"])
    save2output(1, 1, "algae2", algae2["c"])
    save2output(1, 1, "U", variables["U"])
    save2output(1, 1, "Kz", variables["Kz"])
    save2output(1, 1, "C", variables["C"])
    save2output(1, 1, "L", variables["L"])
    save2output(1, 1, "Q2", variables["Q2"])
    save2output(1, 1, "Q2L", variables["Q2L"])
    save2output(1, 1, "N_BV2", variables["N_BV2"])
    save2output(1, 1, "Kq", variables["Kq"])
    save2output(1, 1, "Nu", variables["Nu"])
    #***************************************************************************

    for i in 2:(M-1)

        println("\n\n")

        time = Times[i];
        # println("i = $i")

        # [1] Get forcing variables at time step 
        pressure = get_pressure_at_timestamp(time, Px0, T_Px)
        ustar = calculate_ustar(U)
        W0 = get_wind_speed(i)
        I0 = get_light(i)
        surf = get_surf_flux(i)
        sw = get_shortwave(i) 
        # println("From heat flux ...")
        # println("[HF] Longwave in = ", dfhf[i,"longwave_in"])
        # println("[HF] Longwave out = ", dfhf[i,"longwave_out"])
        println("[HF] Sensible heat = ", dfhf[i,"sensible_heat"])
        println("[HF] Latent heat = ", dfhf[i,"latent_heat"])
        println("[HF] Net longwave radiation = $(dfhf[i,"longwave_in"] - dfhf[i,"longwave_out"])")
        
        # println("Surface water temp= $(C[end])")
        shortwave_in, surface = calculate_heat_flux(forcing_df[i,"DOY"], 
                                                    C[end],
                                                    forcing_df[i,"T_air_C"], 
                                                    forcing_df[i,"shortwave"], 
                                                    forcing_df[i,"RH"], 
                                                    forcing_df[i,"T_dew_point_F"],
                                                    W0)

        println("[HF] shortwave in = ", sw)
        println("[HF] surface flux = ", surf*specific_heat_water*rhoW)
        
        surface_heat_flux = surface/(specific_heat_water * rhoW) 
        
        # [2] Calculate density from temperature field
        rho = calculate_rho(C, base_temp)  
        N_BV2 = calculate_brunt_vaisala(rho, discretization)
        # N_BV2 = clamp.(N_BV2, 0, Inf)

        # [3] Advance velocity field 
        wind_stress = wind_speed_2_wind_stress(W0, discretization) 
        U = advance_velocity(variables, pressure, discretization, wind_stress)

        # [4] Advance TKE / Q2 
        Q2 = advance_Q2(variables, ustar, discretization) 
        Q = @. sqrt(Q2)

        # [5] Advance Q2*L      
        Q2L = advance_Q2L(variables, ustar, discretization)
    
        # [6] Advance temperature 
        body_heat = get_body_flux_heat(shortwave_in, N, z)
        C = advance_scalar(variables, discretization, surface_heat_flux, body_heat) 

        # [7] Semi-implicit: Calculate turbulent lengthscale
        L, Q2L = calculate_lengthscale(Q2, Q2L, N_BV2, discretization)

        # Calculate stability parameter 
        gh = calculate_Gh(N_BV2, L, Q)
        nu_t, Kq, Kz = calculate_turbulent_functions(gh, Q, L, discretization) 

        # [8] Advance phytoplankton
        light = self_shading(algae1, algae2, I0, background_turbidity, discretization)
        gamma = calculate_net_growth(algae1, light, discretization)
        a1 = advance_algae(variables, algae1, gamma, discretization)  # zeros(N) .+ init_algae  #
        algae1["c"] =  clamp.(a1, 1e-5, Inf)   

        # z0_ind = calculate_photic_depth_ind(light, I0)
        # algae["age"] = advance_algae_tracer(variables, algae, gamma, z0_ind, discretization)

        # Algae 2
        # gamma = calculate_net_growth(algae2, light, discretization) 
        # algae2["c"] = advance_algae(variables, algae2, gamma, discretization)
        # algae2["c"] = zeros(N) .+ init_algae 

        # [9] Pack variables for next timestep 
        variables["U"] = U
        variables["C"] = C
        variables["N_BV2"] = N_BV2
        variables["Nu"] = nu_t
        variables["Q2"] = Q2
        variables["Q2L"] = Q2L
        variables["Kq"] = Kq
        variables["Kz"] = Kz
        variables["L"] = L

        if i % isave == 0

            # light_lim = 0.1 * I0
            # search = light .> light_lim
            # ind_photic_depth = length(light) - sum(search) 
            # if ind_photic_depth == 0
            #     ind_photic_depth = 1
            # end
            # println("photic_depth = ", z[ind_photic_depth])

            if i % 5000 == 0
                println("... on time = $i/$M")
            end 
            # println("C= ", C[end-5:end])
            index = div(i, isave) + 1  #(i-1) #div(i, isave)
            save2output(time, index, "algae1", algae1["c"])
            save2output(time, index, "algae2", algae2["c"])
            save2output(time, index, "U", variables["U"])
            save2output(time, index, "Kz", variables["Kz"])
            save2output(time, index, "C", variables["C"])
            save2output(time, index, "L", variables["L"])
            save2output(time, index, "Q2", variables["Q2"])
            save2output(time, index, "Q2L", variables["Q2L"])
            save2output(time, index, "N_BV2", variables["N_BV2"])
            save2output(time, index, "Nu", variables["Nu"])
            save2output(time, index, "Kq", variables["Kq"])
            push!(real_times_saved, real_time[i])
        end

    end

    # ********************** save data ****************************
    units_dict = Dict("U" => "m/s", 
        "C" => "deg C", 
        "Kz" => L"m\$^2\$ s\$^(-1)\$", 
        "algae1" => L"10$^6$/cm$^3$ cells",
        "algae2" => L"10$^6$/cm$^3$ cells",
        "L" => "Turbulent length scale", 
        "Q2" => "TKE", "Q2L" => "TKE*L",
        "N_BV2" => L"s$^(-2)$", "Kq" => "Kq", "Nu" => "Nu_t")

    var2name = Dict("U" => "Velocity", 
                "C" => "Temperature", 
                "Kz" => "Turbulent diffusivity", 
                "algae1" => "Diatom concentration",
                "algae2" => "HAB concentration",
                "L" => "Turbulent length scale", 
                "Q2" => "TKE",
                "Q2L" => "TKE*L",
                "N_BV2" => "Brunt-Vaisala frequency", 
                "Kq" => "Kq", 
                "Nu" => "Turbulent viscosity")


    ds = NCDataset(file_out_name,"c")
    ds.attrib["title"] = "ws1 = $ws1, ws2 = $ws2, pmax1 = $pmax1, pmax2 = $pmax2"

    # model_time = collect(1:M)
    nt = div(M,isave) + 1 
    defDim(ds, "z", length(z)) 
    defDim(ds, "time", nt)

    v = defVar(ds, "z", Float32, ("z",))
    v[:] = z

    v = defVar(ds, "time", Float32, ("time",), attrib = OrderedDict("units" => "seconds"))
    v[:] = collect(1:nt) 

    for var in var2save
        v = defVar(ds, var, Float64,("z","time"), attrib = OrderedDict(
        "units" =>  units_dict[var], "long_name" => var2name[var]))
        v[:,:] = output[var];
    end

    print("Saved $file_out_name \n")
    close(ds)
    
end 


file_out_name = @sprintf("2024_heat_flux.nc") 
run_my_model(1.38e-4, 1.38e-4, 0.008, 0.04, file_out_name)


# using StatProfilerHTML 
# # using ProfileView   
# using Profile 


# @profilehtml run_my_model(ws1, ws2, pmax1, pmax2, file_out_name)

# StatProfilerHTML.view()
# Profile.print() 
