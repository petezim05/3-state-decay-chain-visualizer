include("projFunctions.jl")
using Plots
using Printf

# expects input.txt in format:
# hl_A 1.97
# hl_B 9.91
# Na0 100.0
# Nb0 0.0
# Nc0 0.0
# tfinal 50.0

function read_input(filename)
    params = Dict{String, Float64}()
    open(filename, "r") do f
        for line in eachline(f)
            key, val = split(line)
            params[key] = parse(Float64, val)
        end
    end
    return params
end


#output
function write_output(filename, t, NA, NB, NC, params)
    open(filename, "w") do f
        println(f, "# 3-Component Decay Chain Output")
        println(f, "# Input parameters:")
        println(f, "#   Half-life A (h): $(params["hl_A"])")
        println(f, "#   Half-life B (h): $(params["hl_B"])")
        println(f, "#   Na0: $(params["Na0"])  Nb0: $(params["Nb0"])  Nc0: $(params["Nc0"])")
        println(f, "#   tfinal (h): $(params["tfinal"])")
        println(f, "#")
        println(f, "# t(h)    NA        NB        NC        Total")
        for i in 1:length(t)
            total = NA[i] + NB[i] + NC[i]
            @printf(f, "%.4f   %.6f   %.6f   %.6f   %.6f\n",
                    t[i], NA[i], NB[i], NC[i], total)
        end
    end
end

#so it begins
function main()
    #read
    params = read_input("input.txt")
    λa     = Hl2Dc(params["hl_A"])
    λb     = Hl2Dc(params["hl_B"])
    Na0    = params["Na0"]
    Nb0    = params["Nb0"]
    Nc0    = params["Nc0"]
    tfinal = params["tfinal"]

    #analytical solution
    t_an = collect(0.0:0.01:tfinal)
    NA_an = populationA.(Na0, λa, t_an)
    NB_an = populationB.(Na0, λa, λb, t_an)
    NC_an = populationC.(Na0, λa, λb, t_an)

    #write output
    write_output("output.txt", t_an, NA_an, NB_an, NC_an, params)

    #PLOT 1: NB convergence — coarse, medium, fine + analytical
    dt_values = [1.0, 0.5, 0.25]
    p1 = plot(title="NB Numerical Convergence",
              xlabel="Time (h)", ylabel="Population of B",
              legend=:topright)

    for dt in dt_values
        t_num, N = numAprox(Na0, Nb0, Nc0, λa, λb, dt, tfinal)
        plot!(p1, t_num, N[:,2], label="Numerical Δt = $dt h")
    end
    plot!(p1, t_an, NB_an, label="Analytical", linestyle=:dash, linewidth=2, color=:black)
    savefig(p1, "plot1_convergence.png")

    # PLOT 2: all populations with finest dt
    t_fine, N_fine = numAprox(Na0, Nb0, Nc0, λa, λb, 0.125, tfinal)
    total = N_fine[:,1] .+ N_fine[:,2] .+ N_fine[:,3]

    p2 = plot(t_fine, N_fine[:,1], label="NA",
              title="Decay Chain Populations (Δt = 0.125 h)",
              xlabel="Time (h)", ylabel="Population", legend=:right)
    plot!(p2, t_fine, N_fine[:,2], label="NB")
    plot!(p2, t_fine, N_fine[:,3], label="NC")
    plot!(p2, t_fine, total, label="Total", linestyle=:dash, color=:black)
    savefig(p2, "plot2_populations.png")

    # PLOT 3: t_max NB vs 1/dt, with analytical t_max
    dt_list = [1.0, 0.5, 0.25, 0.125, 0.0625]
    tmax_numerical = Float64[]

    for dt in dt_list
        t_num, N = numAprox(Na0, Nb0, Nc0, λa, λb, dt, tfinal)
        idx = argmax(N[:,2])           # index of maximum NB
        push!(tmax_numerical, t_num[idx])
    end

    tmax_an = t_max_analytical(λa, λb)
    inv_dt = 1.0 ./ dt_list

    p3 = scatter(inv_dt, tmax_numerical,
                 label="Numerical t_max",
                 xlabel="1/Δt (1/h)", ylabel="Time of max NB (h)",
                 title="Convergence of NB Maximum Time",
                 markersize=6)
    hline!(p3, [tmax_an], label="Analytical t_max", linestyle=:dash, color=:red)

    for (x, y, dt) in zip(inv_dt, tmax_numerical, dt_list)
        annotate!(p3, x, y, text("Δt=$(dt)h", :left, 7, :gray30))
    end

    savefig(p3, "plot3_tmax.png")

    println("Done. Output written to output.txt")
    println("Plots saved: plot1_convergence.png, plot2_populations.png, plot3_tmax.png")
end

main()