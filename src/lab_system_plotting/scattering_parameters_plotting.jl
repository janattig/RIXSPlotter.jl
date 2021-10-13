function plot_dq_dependence(ls :: LabSystem, arg, dq_file :: AbstractString, energy_value :: Real, q_beam :: Real=15; plot_label  :: String = "", factor :: Real = 1.0, kwargs...)

    # open file and extract values
    f = open(dq_file, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    vals = [
        get_dq_eps(l, factor = factor) for l in lines
    ]
    close(f)

    # make a list of intensities
    intensities = zeros(length(vals))

    # iterate over all dq values
    for q in 1:length(vals)
        setup_dQ!(ls, vals[q][1], vals[q][2], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls, arg; kwargs...)
        intensities[q] = intensity(spectrum, energy_value)
    end

    # plot the data
    plot(1:length(vals), intensities, label="E = $(energy_value)" * (plot_label=="" ? "" : " "*plot_label))
    xlabel("index")
    ylabel("intensity [a.u.]")
    legend()
end
function plot_dq_dependence(ls :: LabSystem, dq_file :: AbstractString, energy_value :: Real, q_beam :: Real=15; plot_label  :: String = "", factor :: Real = 1.0, kwargs...)

    # open file and extract values
    f = open(dq_file, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    vals = [
        get_dq_eps(l, factor = factor) for l in lines
    ]
    close(f)

    # make a list of intensities
    intensities = zeros(length(vals))

    # iterate over all dq values
    for q in 1:length(vals)
        setup_dQ!(ls, vals[q][1], vals[q][2], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls; kwargs...)
        intensities[q] = intensity(spectrum, energy_value)
    end

    # plot the data
    plot(1:length(vals), intensities, label="E = $(energy_value)" * (plot_label=="" ? "" : " "*plot_label))
    xlabel("index")
    ylabel("intensity [a.u.]")
    legend()
end

function plot_dq_dependence(ls :: LabSystem, arg, dq_file :: AbstractString, energy_values :: Vector{<:Real}, q_beam :: Real=100; plot_label  :: String = "", factor :: Real = 1.0, kwargs...)

    # open file and extract values
    f = open(dq_file, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    vals = [
        get_dq_eps(l, factor = factor) for l in lines
    ]
    close(f)

    # make a list of intensities
    intensities = zeros(length(vals), length(energy_values))

    # iterate over all dq values
    for q in 1:length(vals)
        setup_dQ!(ls, vals[q][1], vals[q][2], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls, arg; kwargs...)
        for i in 1:length(energy_values)
            intensities[q,i] = intensity(spectrum, energy_values[i])
        end
    end

    # plot the data
    for i in 1:length(energy_values)
        plot(1:length(vals), intensities[:,i], label="E = $(energy_values[i])" * (plot_label=="" ? "" : " "*plot_label))
    end
    xlabel("index")
    ylabel("intensity [a.u.]")
    legend()
end
function plot_dq_dependence(ls :: LabSystem, dq_file :: AbstractString, energy_values :: Vector{<:Real}, q_beam :: Real=100; plot_label  :: String = "", factor :: Real = 1.0, kwargs...)

    # open file and extract values
    f = open(dq_file, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    vals = [
        get_dq_eps(l, factor = factor) for l in lines
    ]
    close(f)

    # make a list of intensities
    intensities = zeros(length(vals), length(energy_values))

    # iterate over all dq values
    for q in 1:length(vals)
        setup_dQ!(ls, vals[q][1], vals[q][2], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls; kwargs...)
        for i in 1:length(energy_values)
            intensities[q,i] = intensity(spectrum, energy_values[i])
        end
    end

    # plot the data
    for i in 1:length(energy_values)
        plot(1:length(vals), intensities[:,i], label="E = $(energy_values[i])" * (plot_label=="" ? "" : " "*plot_label))
    end
    xlabel("index")
    ylabel("intensity [a.u.]")
    legend()
end


function plot_dq_dependence(ls :: LabSystem, dq_values, energy_value :: Real, q_beam :: Real=100; plot_label  :: String = "", kwargs...)
    # make a list of intensities
    intensities = zeros(length(dq_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls; kwargs...)
        intensities[q] = intensity(spectrum, energy_value)
    end

    # plot the data
    plot(dq_values ./ pi, intensities, label="E = $(energy_value)" * (plot_label=="" ? "" : " "*plot_label))
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end

function plot_dq_dependence(ls :: LabSystem, arg, dq_values, energy_value :: Real, q_beam :: Real=100; plot_label  :: String = "", kwargs...)
    # make a list of intensities
    intensities = zeros(length(dq_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls, arg; kwargs...)
        intensities[q] = intensity(spectrum, energy_value)
    end

    # plot the data
    plot(dq_values ./ pi, intensities, label="E = $(energy_value)" * (plot_label=="" ? "" : " "*plot_label))
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end

function plot_dq_dependence(ls :: LabSystem, dq_values, i_from :: Int64, i_to :: Int64, q_beam :: Real=100)
    # pass on to other function
    plot_dq_dependence(ls, dq_values, ls.eigensys[:vectors][i_from], ls.eigensys[:vectors][i_to], q_beam, annotate=false)
    # format the plot
    title("intensities for <$(i_to)|D|$(i_from)>")
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end


function plot_dq_dependence(ls :: LabSystem, dq_values, from :: Vector{<:Number}, to :: Vector{<:Number}, q_beam :: Real=100; annotate::Bool=true)

    # print the states
    println("Calculating dq dependence of the following MP states:\n")
    print("Initial "); printMPState(from,basis(ls)); println("")
    print("Final "); printMPState(to,basis(ls)); println("")


    # make a list of intensities
    amplitudes_h  = zeros(Complex{Float64}, length(dq_values))
    amplitudes_v  = zeros(Complex{Float64}, length(dq_values))
    intensities_h = zeros(length(dq_values))
    intensities_v = zeros(length(dq_values))
    intensities   = zeros(length(dq_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        # calculate amplitudes
        amplitudes_h[q]  = get_amplitude(ls.dipole_hor, from,to)
        amplitudes_v[q]  = get_amplitude(ls.dipole_ver, from,to)
        intensities_h[q] = abs.(amplitudes_h[q])^2
        intensities_v[q] = abs.(amplitudes_v[q])^2
        intensities[q]   = intensities_v[q] + intensities_h[q]
    end

    # plot the data
    plot(dq_values ./ pi, real.(amplitudes_h), label="real(A_h)", color=(0,0.2,0.5), alpha=0.4, linestyle="--")
    plot(dq_values ./ pi, imag.(amplitudes_h), label="imag(A_h)", color=(0,0.2,0.5), alpha=0.4, linestyle=":")
    plot(dq_values ./ pi,      intensities_h , label="I_h",       color=(0,0.2,0.5), alpha=1.0, linestyle="-")
    plot(dq_values ./ pi, real.(amplitudes_v), label="real(A_v)", color=(0.7,0,0.1), alpha=0.4, linestyle="--")
    plot(dq_values ./ pi, imag.(amplitudes_v), label="imag(A_v)", color=(0.7,0,0.1), alpha=0.4, linestyle=":")
    plot(dq_values ./ pi,      intensities_v , label="I_v",       color=(0.7,0,0.1), alpha=1.0, linestyle="-")
    plot(dq_values ./ pi,      intensities   , label="I",         color="k"        , alpha=1.0, linestyle="-")
    # format the plot
    if annotate
        title_string = "intensities for <TO|FROM>"
        title(title_string)
        xlabel("dq [pi]")
        ylabel("intensity [a.u.]")
        legend()
    end
end

function plot_dq_dependence(ls :: LabSystem, dq_values, energy_values :: Vector{<:Real}, q_beam :: Real=100; plot_label  :: String = "", kwargs...)
    # make a list of intensities
    intensities = zeros(length(dq_values), length(energy_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls; kwargs...)
        for i in 1:length(energy_values)
            intensities[q,i] = intensity(spectrum, energy_values[i])
        end
    end

    # plot the data
    for i in 1:length(energy_values)
        plot(dq_values ./ pi, intensities[:,i], label="E = $(energy_values[i])" * (plot_label=="" ? "" : " "*plot_label))
    end
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end

function plot_dq_dependence(ls :: LabSystem, arg, dq_values, energy_values :: Vector{<:Real}, q_beam :: Real=100; plot_label  :: String = "", kwargs...)
    # make a list of intensities
    intensities = zeros(length(dq_values), length(energy_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls, arg; kwargs...)
        for i in 1:length(energy_values)
            intensities[q,i] = intensity(spectrum, energy_values[i])
        end
    end

    # plot the data
    for i in 1:length(energy_values)
        plot(dq_values ./ pi, intensities[:,i], label="E = $(energy_values[i])" * (plot_label=="" ? "" : " "*plot_label))
    end
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end
