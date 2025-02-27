# Library to deal with model output 
output = Dict()
times = [] 

# ***********************************************************************
function create_output_dict(M::Int, isave::Int, Vars::Vector{<:String}, N::Int, output=output)
    nvars = length(Vars)
    n_saved_steps = div(M,isave) + 1 #- 1

    for i = 1:nvars
        output[Vars[i]] = zeros(Float64, N, n_saved_steps)
    end
end 
# ***********************************************************************


# ***********************************************************************
function save2output(time, index, varname, value, output=output)
    output[varname][:,index] .= value
    push!(times, time)
end 
# ***********************************************************************

function plot_variable(varname, discretization, output=output)
    ntimes = length(unique(times))
    print(output)
    print("Plotting $(ntimes) steps")
    z = discretization["z"]
    p = plot() 

    for i in 1:(ntimes-1)
        print(i)
        print(output[varname][:,i])
        plot!( output[varname][:,i], -z) # label="t = $(times[i])")
    end
   display(p) # plot(times, output[varname][1,:])
end