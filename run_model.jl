#!/usr/bin/env julia

using Plots
using Printf
using DataStructures: OrderedDict
using NCDatasets
using Arrow, DataFrames
using CSV, DataFrames
# using Colors
# using ColorSchemes
# using Plots
using Printf
using LaTeXStrings
using Profile
using Statistics 

include("model_code/calculate_physical_variables.jl") 
include("model_code/advance_variables.jl")
include("model_code/phytoplankton.jl")
include("model_code/forcings.jl") 
include("model_code/output.jl")
include("model_code/calculate_heat_fluxes.jl")
include("model_code/kaimal.jl")

# ws1 = parse(Float64, ARGS[1])
# ws2 = parse(Float64, ARGS[2])
# pmax1 = parse(Float64, ARGS[3])
# pmax2 = parse(Float64, ARGS[4])
# fout_name = ARGS[5]



function run_my_model(Nens::Int, file_out_name::String)

    @info "Running model with $Nens ensemble members and saving to $file_out_name ..."
    forcing_fn = "/global/homes/s/siennaw/scratch/siennaw/stockton_field_data/forcing_for_model/2024/august6-28/"

    # CIMIS Data 
    forcing_fp = @sprintf("%s/CIMIS2.csv", forcing_fn)  
    @info ("Reading in forcing data from $forcing_fp ...")
    forcing_df = CSV.read(forcing_fp, DataFrame)

    # Heat flux Data 
    heat_flux_fn =  @sprintf("%s/heat_flux2.csv", forcing_fn)  
    @info ("Reading in heat flux data from $heat_flux_fn ...")
    dfhf = CSV.read(heat_flux_fn, DataFrame)


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



    #********************** SPATIAL DOMAIN  ***************************
    N = 60    # number of grid points
    H = 6    # depth (meters)
    dz = H/N  # grid spacing - may need to adjust to reduce oscillations
    dt = 10   # (seconds) size of time step 
    M  = 360*24 #190050 #359*24*22 #206000 # 359*24* 25 # 17280 #00 #000 # 50000  #500 # 17280 # 
            

    # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
    isave = 100 # 360 #1000
    var2save = ["U","Kq", "Nu", "C", "Kz", "L", "Q2", "Q2L", "N_BV2", "algae1","algae2"]

    create_output_dict(M, Nens, isave, var2save, N)

    # Create depth vector 
    z = collect(H:-dz:dz) .- dz/2 # depth vector


        #***********************************************************************
    # Wind time series 
    # wind_fn =  @sprintf("%s/wind.csv", forcing_fn)  
    # @info ("Reading in wind data from $wind_fn ...")
    # df = CSV.read(wind_fn, DataFrame)
    # wind = df[!,"WindSpeed"] 
    # real_time = df[!,"time"]
    Uwind = 3.5 


        
        # function get_wind_speed(index::Int, wind=wind)
        #     return wind[index]
        # end
        # #***********************************************************************


        
    function get_body_flux_heat(sw, N, z)
        rad = zeros(N)       # radiation 
        qsource = zeros(N)   # heat source term
        # A = 0.78             # From Paulson & Simpson for Type III water 
        # g1 = 1.4
        # g2 = 7.9

        # new coeffiicents 3/13 (type 9 water)
        g1 = 0.25
        g2 = 0.55
        A = 0.91

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


        # (4) Light 
        DIURNAL_LIGHT = false  
        background_turbidity =  0.5

        #********************** DEFINE PHYTOPLANKTON FORCINGS ***************************
        # init_algae = 0.005
        init_algae = 1

        algae1 = Dict("k" => 0.034,              # specific light attenuation coefficient [cm^2 / 10^6 cells]
        "pmax" =>0,           # maximum specific growth rate [1/hour]
        "ws" => 1.38e-4, #1.38e-4,           # vertical velocity [m/s]
        "Hi" => 40,                # half-saturation of light-limited growth [mu mol photons * m^2/s]
        "Li" => 0.004 * hr2s,             # specific loss rate [1/hour]
        "name" => "HAB",           # name of the species
        "self_shading" => true)    # self-shading effect (true/false)

        algae2 = Dict("k" => 0.034,              # specific light attenuation coefficient [cm^2 / 10^6 cells]
                    "pmax" => 0, #0.008 * hr2s,           # maximum specific growth rate [1/hour]
                    "ws" => -1.38e-4, #1.38e-4,           # vertical velocity [m/s]
                    "Hi" => 40,                # half-saturation of light-limited growth [mu mol photons * m^2/s]
                    "Li" => 0.004 * hr2s,             # specific loss rate [1/hour]
                    "name" => "HAB",           # name of the species
                    "self_shading" => true)    # self-shading effect (true/false)

        #***************************************************************************
        #   Initialize variables
        #***************************************************************************

        # Create dictionary to hold important discretization parameters
        discretization = Dict("beta" => (dt/dz^2), "dz" => dz, "dt" => dt, "N" => N, "z"=> z, "H" => H)
        wind_out = zeros(Nens, div(M,isave) + 1)

        for Ne in 1:Nens

            println("Running ensemble member $Ne / $Nens ...")

            wt, wind = simulate_kaimal_wind(M, 1/dt, Uwind)
            wind_out[Ne, :] = wind[1:isave:end]
            # Initalize velocity
            U = similar(z) .+ 1e-1

            # Initalize temperature field
            C = [ 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.786400, 27.798220, 27.816607, 27.834993, 27.853380, 27.878350, 27.915170, 27.951990, 27.988810, 28.025630, 28.062450, 28.099270, 28.136090, 28.172910, 28.181460, 28.182300, 28.183140, 28.183980, 28.184820, 28.185660, 28.186500, 28.187340, 28.188050, 28.188283, 28.188517, 28.188750, 28.188983, 28.189087, 28.189180, 28.189273, 28.189367, 28.189280, 28.189093, 28.188907, 28.188720, 28.188400, 28.187840, 28.187280, 28.186720, 28.186233, 28.186700, 28.187167, 28.187027, 28.186840, 28.188853, 28.191467, 28.189520, 28.185040,]

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
            save2output(1, 1, Ne, "algae1", algae1["c"])
            save2output(1, 1, Ne, "algae2", algae2["c"])
            save2output(1, 1, Ne, "U", variables["U"])
            save2output(1, 1, Ne, "Kz", variables["Kz"])
            save2output(1, 1, Ne, "C", variables["C"])
            save2output(1, 1, Ne, "L", variables["L"])
            save2output(1, 1, Ne, "Q2", variables["Q2"])
            save2output(1, 1, Ne, "Q2L", variables["Q2L"])
            save2output(1, 1, Ne, "N_BV2", variables["N_BV2"])
            save2output(1, 1, Ne, "Kq", variables["Kq"])
            save2output(1, 1, Ne, "Nu", variables["Nu"])
            #***************************************************************************

            for i in 2:(M-1)

                time = Times[i];


                # [1] Get forcing variables at time step 
                pressure = get_pressure_at_timestamp(time, Px0, T_Px)
                ustar = calculate_ustar(U)
                W0 = wind[i] #get_wind_speed(i) 
                W0 = clamp(W0, 0, Inf)
                I0 = get_light(i)
                surf = get_surf_flux(i)
                sw = get_shortwave(i) 

                shortwave_in, surface = calculate_heat_flux(forcing_df[i,"DOY"], 
                                                            C[end],
                                                            forcing_df[i,"T_air_C"], 
                                                            forcing_df[i,"shortwave"], 
                                                            forcing_df[i,"RH"], 
                                                            forcing_df[i,"T_dew_point_F"],
                                                            W0)

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
                # light = self_shading(algae1, algae2, I0, background_turbidity, discretization)
                gamma = zeros(N) #calculate_net_growth(algae1, light, discretization)
                a1 = advance_algae(variables, algae1, gamma, discretization)  # zeros(N) .+ init_algae  #
                algae1["c"] =  clamp.(a1, 1e-5, Inf)   
                # Algae 2
                # gamma = calculate_net_growth(algae2, light, discretization) 
                a2 = advance_algae(variables, algae2, gamma, discretization)
                algae2["c"] = clamp.(a2, 1e-5, Inf)   

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

                    if i % 5000 == 0
                        println("... on time = $i/$M")
                    end 
                    # println("C= ", C[end-5:end])
                    index = div(i, isave) + 1  #(i-1) #div(i, isave)
                    save2output(time, index, Ne, "algae1", algae1["c"])
                    save2output(time, index, Ne, "algae2", algae2["c"])
                    save2output(time, index, Ne, "U", variables["U"])
                    save2output(time, index, Ne, "Kz", variables["Kz"])
                    save2output(time, index, Ne,  "C", variables["C"])
                    save2output(time, index, Ne, "L", variables["L"])
                    save2output(time, index, Ne, "Q2", variables["Q2"])
                    save2output(time, index, Ne, "Q2L", variables["Q2L"])
                    save2output(time, index, Ne, "N_BV2", variables["N_BV2"])
                    save2output(time, index, Ne, "Nu", variables["Nu"])
                    save2output(time, index, Ne, "Kq", variables["Kq"])
                end
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
        # ds.attrib["title"] = "ws1 = $ws1, ws2 = $ws2, pmax1 = $pmax1, pmax2 = $pmax2"

        # model_time = collect(1:M)
        nt = div(M,isave) + 1 
        defDim(ds, "z", length(z)) 
        defDim(ds, "time", nt)
        defDim(ds, "Nens", Nens)

        v = defVar(ds, "z", Float32, ("z",))
        v[:] = z

        v = defVar(ds, "time", Float32, ("time",), attrib = OrderedDict("units" => "seconds"))
        v[:] = collect(1:nt) 

        v = defVar(ds, "wind" , Float64,("Nens", "time",), attrib = OrderedDict(
                "units" =>  "m/s", "long_name" => "generated wind speed"))
        v[:] = wind_out #wind[1:isave:end];

        for var in var2save
            v = defVar(ds, var, Float64,("Nens","z","time"), attrib = OrderedDict(
            "units" =>  units_dict[var], "long_name" => var2name[var]))
            v[:,:,:] = output[var];
        end

        print("Saved $file_out_name \n")
        close(ds)
     
        
end 


file_out_name = @sprintf("2024_heat_flux_ens100.nc") 
run_my_model(100, file_out_name)


    # using StatProfilerHTML 
    # # using ProfileView   
    # using Profile 


    # @profilehtml run_my_model(ws1, ws2, pmax1, pmax2, file_out_name)

    # StatProfilerHTML.view()
    # Profile.print() 
