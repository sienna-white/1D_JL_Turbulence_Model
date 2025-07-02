# CONSTANTS 
altitude = 4    # Altitude of Stockton
latitude = 37.5 # Latitude in degrees (e.g., 37.5 for California)

# [1] for shortwave
Tropic = 23.45 
days_per_year = 365.0       
DclDay = 284.0
DgCrcl = 360.0 
deg2rad = 57.29577951 
RefInd = 1.33                        # This is a mild function of T and S. http//scubageek.com/articles/wwwh2o.html
latitude_radian = latitude/deg2rad

# [2] for longwave
surface_pressure = 101325 * (1 - altitude * 2.25577e-5)^5.25588 * (1/100) # calculated in (mb)
surface_pressure = surface_pressure * 100 # Convert to Pascal


function calculate_heat_flux(DOY, T_water_C, T_air_C, shortwave, RH, T_dew_point_F, wind_speed)
    println("incoming shortwave radiation (W/m^2) = $(shortwave)")
    # [1] Saturation vapor pressure (hPa)
    es = get_saturation_vapor_pressure(T_air_C) 
    es = es * 100 # Conversion millibar --> Pascal

    # [2] Actual water vapor pressure (hPa)
    ea = RH/100 * es

    # [3] Calculate zenith 
    zenith = calculate_zenith(DOY)

    shortwave_in = calculate_shortwave_in(DOY, zenith, shortwave)  # Shortwave radiation in (W/m^2)
    sensible_heat_flux, latent_heat_flux  = get_sensible_and_latent(wind_speed, T_water_C, T_air_C, es, ea)  # Sensible and latent heat fluxes (W/m^2)
    net_longwave = calculate_net_longwave(T_water_C, T_air_C, ea)

    println("T_water_C = $(T_water_C)")
    println("T_air_C = $(T_air_C)")
    println("Sensible heat flux = $(sensible_heat_flux)")
    println("Latent heat flux = $(latent_heat_flux)")
    surface = net_longwave - sensible_heat_flux - latent_heat_flux  # Surface energy balance
    println("Shortwave radiation = $(shortwave_in)")
    println("Surface energy balance = $(surface)")   

    return shortwave_in, surface
end 
function calculate_shortwave_in(DOY, zenith, shortwave)
    if shortwave > 10
        albedo = calculate_albedo(zenith)
        Qsin = shortwave * (1 - albedo)
    else 
        Qsin = 0.0  # If shortwave is less than 5 W/m^2, set it to zero
    end
    return Qsin
end 

function calculate_net_longwave(T_water_C, T_air_C, ea)
    ccf = 0.665788  # Cloud cover fraction (from the original code)
    cloud = 0.1
    epsilon_w = 0.972      # Emissivity of water (Davies et al., 1971)
    sigma    = 5.67e-8     # Stefan-Boltzmann constant (W m^-2 K^-4)
    x1 = (1.0-ccf*cloud*cloud)*(celsius2kelvin(T_air_C)^4)
    x2 = (0.39 - 0.05*sqrt(ea*0.01))
    x3=4.0*(celsius2kelvin(T_air_C)^3) * (T_water_C - T_air_C)
    ql = -epsilon_w*sigma*(x1*x2 + x3)
    return ql 
end 


function celsius2kelvin(celsius)
    """Convert Celsius to Kelvin."""
    return celsius + 273.13
end 


function calculate_zenith(DOY, DclDay=DclDay, DgCrcl=DgCrcl, Tropic=Tropic, deg2rad=deg2rad, latitude=latitude)
    deg2rad = pi / 180
    rad2deg = 180 / pi

    # 37.9283,-121.2893
    latitude_r  = 37.9283 * deg2rad  # Convert latitude to radians
    longitude_r = -121.2893 * deg2rad  # Convert longitude to radians
    
    days_per_year = 365.25
    th0 = 2*pi*DOY/days_per_year
    th02 = 2*th0
    th03 = 3*th0
    
    # sun declination :
    sundec = (0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)       
            - 0.006758*cos(th02) + 0.000907*sin(th02)                 
            - 0.002697*cos(th03) + 0.001480*sin(th03)) 

    # sun hour angle 
    hour  = (DOY-floor(DOY))*24
    thsun = (hour-12)*15*deg2rad + longitude_r   

    # !  cosine of the solar zenith angle :
    coszen =sin(latitude_r)*sin(sundec) + cos(latitude_r)*cos(sundec)*cos(thsun)
    if (coszen < 0) 
        coszen = 0
    end 

    solar_zenith_angle = acos(coszen)
    return solar_zenith_angle 
end 


function calculate_albedo(zenith)
    zenith = clamp(zenith, 0.0, 1.48353) # clamp around pi/2 (85 deg)
    RefInd = 1.33      
    # Angle of Refraction calculation based on Snell's Law 
    RefAng = asin(sin(zenith)/RefInd)  # Angle of refraction
    albedo1 = 0.026 + 0.065 * (tan(zenith))^1.8
    println("Albedo1 = $albedo1")
    
    # Albedo Calculation 
    A1 = tan(zenith - RefAng)^2
    A2 = tan(zenith + RefAng)^2
    A3 = sin(zenith - RefAng)^2
    A4 = sin(zenith + RefAng)^2
    albedo = 0.5 * (A1/A2 + A3/A4)
    println("Albedo2 = $albedo")
    albedo = clamp(albedo, 0.02, 0.2) # 0.5 worked well 
    return albedo 
end 

function get_saturation_vapor_pressure(T_air_C)
    # Empirical fit from P.R. Lowe, 1976. An Approximating Polynomial for the
    # Computation of Saturation Vapor Pressure. Journal of Applied Meteorology
    a0 = 6.107799961
    a1 = 4.436518521e-1
    a2 = 1.428945805e-2
    a3 = 2.650648471e-4
    a4 = 3.031240396e-6
    a5 = 2.034080948e-8
    a6 = 6.136820929e-11

    # Saturation vapor pressure of the air 
    sat_pressure = a0 + T_air_C *(a1 + T_air_C*(a2 + T_air_C*(a3 + T_air_C*(a4 + T_air_C*(a5 + a6*T_air_C)))))
    return sat_pressure
end 


function get_kondo_coeff2(wind_speed)
    ae_h = ( 0., 0.927, 1.15 , 1.17 , 1.652 )
    ae_e = ( 0., 0.969, 1.18 , 1.196, 1.68  )

    be_h = ( 1.185, 0.0546, 0.01  , 0.0075, -0.017 )
    be_e = ( 1.23 , 0.0521, 0.01  , 0.008 , -0.016 )

    ce_h = ( 0., 0., 0., -0.00045, 0. )
    ce_e = ( 0., 0., 0., -0.0004 , 0. )

    pe_h = ( -0.157, 1., 1., 1., 1. )
    pe_e = ( -0.16 , 1., 1., 1., 1. )

    eps = 1.0e-12
    if wind_speed < 2.2
        x = log(wind_speed + eps)
        chd = (be_h[1]*exp(pe_h[1]*x))*1.0e-3
        ced = (be_e[1]*exp(pe_e[1]*x))*1.0e-3
    elseif wind_speed < 5
        x = exp(log(wind_speed+eps))
        chd = (ae_h[2]+be_h[2]*x)*1.0e-3
        ced = (ae_e[2]+be_e[2]*x)*1.0e-3
    elseif wind_speed < 8
        x = exp(log(wind_speed+eps))
        chd = (ae_h[3]+be_h[3]*x)*1.0e-3
        ced = (ae_e[3]+be_e[3]*x)*1.0e-3
    elseif wind_speed < 25
        x = exp(log(wind_speed+eps))
        chd = (ae_h[4]+be_h[4]*x+ce_h[4]*(wind_speed-8.0)^2)*1.0e-3
        ced = (ae_e[4]+be_e[4]*x+ce_e[4]*(wind_speed-8.0)^2)*1.0e-3
    end 
    return chd, ced
end 

function get_sensible_and_latent(wind_speed, T_water_C, T_air_C, es, ea, surface_pressure=surface_pressure)
    SpecificHeatAir = 1006.43 # Specific heat capacity of air, J kg-1 K-1
    GasConstant = 287.1   # gas constant for dry air J kg-1 K-1
    HumidityRatio=0.62198        # Constant for specific humidity calculation
    eps = 1.0e-12

    # Latent heat of vaporization (J kg^-1)
    Lv = (2.5 - 0.00234*T_water_C)*1e6
    
    # Specfic humidity (kg kg^-1)
    qa = HumidityRatio*ea/(surface_pressure-0.377*ea)

    # Saturation specific humidity
    qs = HumidityRatio*es/(surface_pressure-0.377*es) 

    # Density of the air  
    rho_air = surface_pressure/(GasConstant*celsius2kelvin(T_air_C)*(1.0 + HumidityRatio*qa))
    
    s0 = 0.25*(T_water_C - T_air_C)/(wind_speed+1e-10)^2
    s = s0*abs(s0)/(abs(s0)+0.01)

    chd, ced = get_kondo_coeff2(wind_speed)  # Get Kondo coefficients based on wind speed
    # println("[2nd method] Kondo coefficients: chd = $chd, ced = $ced")
    # Adjust coefficients based on stability
    # println("Stability (s) = $s")
    if s < 0
        x = 0 
        if s > -3.3
            x = 0.1+0.03*s + 0.9*exp(4.8*s)
        else
            x = 0 # 
        end 
        chd = chd * x
        ced = ced * x

    else
        chd=chd*(1.0 + 0.63 * sqrt(s))
        ced=ced*(1.0 + 0.63 * sqrt(s))
    end 
    # Latent heat flux
    ced = clamp(ced, 0.0022, 1)
    latent = ced * Lv * rho_air * wind_speed * (qs - qa) 

    # Sensible heat flux 
    sensible = ced * SpecificHeatAir * rho_air * wind_speed *(T_water_C - T_air_C)      
    # println("Wind speed = $wind_speed m/s")
    # println("ced = $ced, Lv = $Lv, rho_air=$rho_air, Δq = $(qs - qa)")
    # println("Final ced = $ced")

    return sensible, latent
end 



    # Kondo coefficients for different wind speed ranges
    # ad, ah, ae, bd, bh, be, ch, ce, pd, ph, pe = get_kondo_coeff(wind_speed)
    # cdd = (ad + bd * exp(pd * log(wind_speed + eps)))*1e-3
    # chd = (ah + bh * exp(ph * log(wind_speed + eps)) + ch * (wind_speed-8.0)^2)*1e-3
    # ced = (ae + be * exp(pe * log(wind_speed + eps)) + ce * (wind_speed-8.0)^2)*1e-3
    # println("[1st method] Kondo coefficients: chd = $chd, ced = $ced")



#     sat_pressure = get_saturation_vapor_pressure(T_air_C)   # Saturation vapor pressure in hPa
# vapor_pressure = RH/100 * sat_pressure  
    # Actual water vapor pressure
    # ea = get_saturation_vapor_pressure(T_air_C)
    # [2] Specific humidity at saturation pressure (kg kg^-1)  
    # q0 = 0.622 * e_sat /surface_pressure      
    # specific_humidity = HumidityRatio*ea/(surface_pressure-0.377*ea)
    
    

# sat_pressure = get_saturation_vapor_pressure(T_air_C)   # Saturation vapor pressure in hPa
# vapor_pressure = RH/100 * sat_pressure                  # Vapor pressure 
    
    # print("specific_humidity = %2.2f kg/kg" % specific_humidity)

    # specific_humidity = 0.622 * vapor_pressure / surface_pressure  # Specific humidity of the air at height z_q above the water surface (kg kg^-1)
    # print("specific_humidity = %2.2f kg/kg" % specific_humidity)

    # humidity_at_saturation = 0.622 * e_sat / surface_pressure 

    # Calculate density of the air 
    # Ra = GasConstant * (1 + 0.608*specific_humidity)


    # rho_air = 100 * surface_pressure/(Ra*(T_air_C + 275.16))
    # print("density of air (rho_air) = %2.2f kg/m^3" % rho_air)

        # [1] Saturated vapor pressure (hPa) 
    # e_sat = 6.11 * exp(17.27*T_water_C/(237.3 + T_water_C)) 

    # ea = get_saturation_vapor_pressure(T_air_C)   # Saturation vapor pressure in hPa



# function calculate_dew_point(vapor_pressure)
#     # """Calculate the dew point temperature from vapor pressure."""
#     Td1 = 243.5 * log(vapor_pressure / 6.112)
#     Td2 = 17.67 - log(vapor_pressure / 6.112)
#     T_dew_point = Td1 / Td2 # Degrees Celsius 
#     T_dew_point_F = T_dew_point*9/5 + 32 # Convert to Fahrenheit
#     return T_dew_point_F

# T_air_k =  celsius2kelvin(T_air_C)  # air temperature in Kelvin
# T_water_k  = celsius2kelvin(T_water_C)  # Convert surface temperature to Kelvin
# sat_pressure = get_saturation_vapor_pressure(T_air_C)   # Saturation vapor pressure in hPa
# vapor_pressure = RH/100 * sat_pressure                  # Vapor pressure 
# T_dew_point_F = calculate_dew_point(vapor_pressure)  # Dew point temperature in Fahrenheit
# print("The dew point temperature is #2.2f F" % T_dew_point_F)


# function calculate_longwave_out(T_water_C)
#     # Calculate longwave radiation out (W/m^2) using the Stefan-Boltzmann law
#     T_water_k = celsius2kelvin(T_water_C)  # Convert surface temperature to Kelvin
#     epsilon_w = 0.972      # Emissivity of water (Davies et al., 1971)
#     sigma    = 5.67e-8     # Stefan-Boltzmann constant (W m^-2 K^-4)
#     Qlout = epsilon_w * sigma * T_water_k^4   
#     return Qlout 
# end 


# function calculate_longwave_in(DOY, shortwave, T_dew_point_F, T_air_C, ea, zenith, surface_pressure=surface_pressure) 
#     G = 2.92     # Values from Smith-Gamma table (need to look up)# Other option for springtime is 3.11 
#     cos_zenith = cos(zenith)                             # Cosine of the solar zenith angle
#     air_mass_thickness_coeff = 35 * (1244 * cos_zenith^2 + 1)^(-1/2) 
#     precip_water = exp((0.1133 - log(G+1)) + 0.0393*T_dew_point_F) # Note T_dew_point is in Fahrenheit 

#     # Aerosol attenuation calculated from clear-sky conditions 
#     T_a = 0.95^air_mass_thickness_coeff   # Houghton (1954)

#     Tr_x_Tpg = 1.021 - 0.084*sqrt(air_mass_thickness_coeff*(0.000949*surface_pressure + 0.0151))

#     # Water vapor absorption in the atmosphere from McDonald (1960) 
#     T_w = 1 - 0.077*(precip_water * air_mass_thickness_coeff)^(0.3)

#     I_eff = 1353 * (1 + 0.034 * cos(2*pi * (DOY-1)/365))  # (Meyers and Dale, 1983) # Sw note cos might be^2 ? 

#     # Calculate cloud cover fraction 
#     clear_sky_shortwave =  I_eff * cos_zenith * Tr_x_Tpg * T_w * T_a # Clear-sky shortwave radiation (W m^-2)

#     s = shortwave/clear_sky_shortwave  # Ratio of measured shortwave radiation to the clear-sky shortwave radiation
#     clf = clamp((1 - s), 0.0, 0.5)  # Cloud cover fraction 
#     println("Cloud cover fraction (clf) = $clf (compare to 0.665788 )")

#     sigma    = 5.67e-8   
#     month = 8 
#     term1 = sigma * celsius2kelvin(T_air_C)^4 
#     term2 = 1.22 + 0.06 * sin((month + 2)* pi/6)
#     term3 = (ea/celsius2kelvin(T_air_C))^(1/7)
#     emissivity = (clf + (1-clf) * term2 * term3)
#     emissivity = clamp(emissivity, 0.0, 1)  # Ensure emissivity is between 0 and 1
#     Qlin = term1 * emissivity #(clf + (1-clf) * term2 * term3) 
#     # println("term1 = $term1, term2 = $term2, term3 = $term3, ea = $ea, T_air_C = $T_air_C")

#     return Qlin 
# end 




# [3] for heat flux
# C_pa = 1006 # Specific heat of air at constant pressure (J kg^-1 C^-1)

    # longwave_out = calculate_longwave_out(T_water_C)  # Longwave radiation out
    # longwave_in  = calculate_longwave_in(DOY, shortwave, T_dew_point_F, T_air_C, ea, zenith)  # Longwave radiation in (W/m^2)


        # z = DclDay + floor(DOY)
    # x = DgCrcl * z/days_per_year/deg2rad
    # y = sin(x)
    # Decl = Tropic/deg2rad * y

    # # Hour-angle calculation where time is hours from midnight  
    # HrAng = ((DOY-floor(DOY))*24 -12)*15.0 * deg2rad # change divide to multiply 
 
    # # Zenith angle calculation 
    # zenith = acos(sin(Decl)*sin(latitude_radian)+cos(Decl)*cos(latitude_radian)*cos(HrAng))
    # println("Zenith angle (radians) = $zenith")

        # Calculate hortwave radiation in (W/m^2)
    # zenith = calculate_zenith(DOY)
        # albedo = calculate_albedo(zenith)
        # println("Albedo = $albedo")
        # Qsin = shortwave * (1 - albedo)



    # φ = 37.5 * deg2rad
    
    # # # Approximate solar declination (Cooper 1969)
    # δ = 23.44 * deg2rad * sin(2*pi * (DOY - 81) / 365)

    # # # Solar hour angle (radians)
    # h = ((DOY-floor(DOY))*24 - 12.0) * 15.0 * deg2rad

    # # # Solar zenith angle
    # cos_zenith = sin(φ)*sin(δ) + cos(φ)*cos(δ)*cos(h)
    # # cos_zenith = clamp(cos_zenith, -1.0, 1.0)  # avoid domain errors
    # zenith = acos(cos_zenith)  # in radians

# function get_kondo_coeff(wind_speed)
#     if wind_speed < 2.2
#         ad = 0
#         ah = 0 
#         ae = 0
#         bd = 1.08
#         bh = 1.185
#         be = 1.23
#         ch = 0 
#         ce = 0 
#         pd = -0.15
#         ph = -0.157
#         pe = -0.16
#     elseif wind_speed < 5 
#         ad = 0.771
#         ah = 0.927
#         ae = 0.969
#         bd = 0.0858
#         bh = 0.0546
#         be = 0.0521
#         ch = 0 
#         ce = 0 
#         pd = 1
#         ph = 1
#         pe = 1
#     elseif wind_speed < 8   
#         ad = 0.867
#         ah = 1.15
#         ae = 1.18
#         bd = 0.0667
#         bh = 0.01
#         be = 0.01
#         ch = 0 
#         ce = 0 
#         pd = 1
#         ph = 1
#         pe = 1
#     elseif wind_speed < 25   
#         ad = 1.2
#         ah = 1.17
#         ae = 1.196
#         bd = 0.025
#         bh = 0.0075
#         be = 0.008
#         ch = -0.00045
#         ce = -0.0004
#         pd = 1
#         ph = 1
#         pe = 1
#     end 
#     return ad, ah, ae, bd, bh, be, ch, ce, pd, ph, pe
# end 

