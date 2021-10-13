rc("figure", dpi=200)
rc("font",  family="serif")
rc("xtick", labelsize="x-small")
rc("xtick", direction="in")
rc("ytick", labelsize="x-small")
rc("ytick", direction="in")
rc("lines", linewidth=1)
rc("figure", dpi=200)



################################################################################
#
#   Spectrum plotting functions
#
################################################################################


# plot a spectrum that has been calculated before
function plot_spectrum(
        spectrum    :: S,
        energies    :: Vector{<:Real};
        new_figure  :: Bool = true,
        show_figure :: Bool = true,
        plot_label  :: String = "",
        plot_color  :: Any = "b",
        plot_transitions :: Bool = false
    ) where {T,S<:AbstractSpectrum{T}}

    # configure the plot
    if new_figure
        figure()
        xlabel("energy [meV]")
        ylabel("intensity [a.u.]")
    end

    # plot the spectrum
    plot(
        energies, [intensity(spectrum, omega) for omega in energies],
        label=plot_label, color=plot_color
    )

    # add a legend
    if plot_label != ""
        legend()
    end

    # maybe plot transitions
    if plot_transitions
        transitions = spectrum.transitions
        intensities = weight.(transitions)
        intensities ./= sum(intensities)
        intensities ./= maximum(intensities)
        for t in 1:length(transitions)
            if intensities[t] > 0.01
                axvline(frequency(transitions[t]), color="k", alpha=intensities[t])
            end
        end
    end

    # set the x limits
    xlim(energies[1], energies[end])
    ylim(0, ylim()[2])

    # tighten the layout
    tight_layout()

    # show the plot
    if show_figure
        show()
    end
end




# export functions
export plot_spectrum