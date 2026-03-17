using FFTW
using Random


"""
Calculates the Kaimal Power Spectral Density.
"""
function kaimal_spectrum(f, z, U, u_star)
    # Avoid division by zero by replacing 0 with a small epsilon
    f_safe = ifelse.(f .== 0, 1e-6, f)
    
    n = (f_safe .* z) ./ U
    S = (u_star^2 * 105 .* n) ./ (f_safe .* (1 .+ 33 .* n).^(5/3))
    return S
end

"""
Simulates wind time series.
"""
function simulate_kaimal_wind(M, fs, U)

    # """
    # Simulates wind time series.
    # T: Total time (s)
    # fs: Sampling frequency (Hz)
    # z: Height (m)
    # U: Mean wind speed (m/s)
    # I: Turbulence intensity (decimal, e.g., 0.15 for 15%)

    
    z = 30          # meters
    # mean_U = 8         # m/s
    intensity = 1    # 15% turbulence
    I = intensity 


    N = M #int(T * fs)
    dt = fs / N

    # @info "Generating Kaimal wind time series with T=$M s, fs=$fs Hz, z=$z m, U=$U m/s, I=$(I*100)% turbulence intensity..."
    
    # N = floor(Int, T * fs)
    # dt = 1/fs
    
    # Frequency vector for rFFT (0 to Nyquist)
    # Julia's rfftfreq equivalent:
    f = collect(range(0, stop=fs/2, length=floor(Int, N/2) + 1))
    
    # Turbulence intensity & friction velocity
    sigma_u = I * U
    u_star = sigma_u / 2.45 
    
    # Calculate target PSD
    S_target = kaimal_spectrum(f, z, U, u_star)
    
    # Generate random phases
    # exp(im * theta) is the Julia equivalent of exp(1j * theta)
    random_phases = exp.(im .* 2π .* rand(length(f)))
    
    # Construct complex spectrum
    # Using the same scaling for irfft as the Python snippet
    Z = sqrt.(S_target .* fs ./ 2) .* random_phases
    
    # Inverse FFT to get time series
    # irfft in Julia requires the length of the output signal if N is odd
    u_prime = irfft(Z, N)
    
    # Add mean wind speed
    u = U .+ u_prime.*20
    
    t=0
    # t = range(0, M, length=N)
    u[1:(360*14)] .= u[1:(360*14)] .* 0.5
    return t, u
end

# # --- Usage ---
# # Parameters
# T_total = 600       # 10 minutes
# fs_rate = 1/10      # output frequency (in Hz)
# height = 2          # meters
# mean_U = 8          # m/s
# intensity = 0.2     # 20% turbulence

# t, u = simulate_kaimal_wind(M, T, fs, z, U, I)

# plot(t, u, title="Simulated Kaimal Wind Series", xlabel="Time (s)", ylabel="Wind Speed (m/s)", legend=false)