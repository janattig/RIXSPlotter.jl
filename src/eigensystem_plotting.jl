function show_operator(op :: AbstractOperator)
    m = matrix_representation(op)
    figure(figsize=(9,4))
    subplot(121)
    title("real part")
    imshow(real.(m)/maximum(real.(m)), cmap="RdBu", vmin=-1.1, vmax=1.1)
    subplot(122)
    title("imaginary part")
    imshow(imag.(m)/maximum(imag.(m)), cmap="RdBu", vmin=-1.1, vmax=1.1)
end
export show_operator


function show_energy_evolution(op :: AbstractOperator, parameter::Symbol, values; site=:all, subtract_GS::Bool = false, figsize::Tuple=(8,5), new_figure=true, dumpfile::String="", color="b", kwargs...)
    if new_figure
        figure(figsize=figsize)
    end
    bands = [zeros(length(values)) for i in 1:length(basis(op))]
    ba_param = get_parameter(op, parameter, site=site)
    for i in 1:length(values)
        set_parameter!(op, parameter, values[i]; recalculate=true, site=site)
        evals = energies(op)
        for j in 1:length(evals)
            bands[j][i] = evals[j]
        end
        if subtract_GS
            for j in 1:length(evals)
                bands[j][i] -= evals[1]
            end
        end
    end
    # saving
    if dumpfile != ""
        # open file
        f = open(dumpfile, "w")
        # write header line
        lines = split(string(op), "\n")
        for l in lines
            print(f,"# ",l, "\n")
        end
        print(f, "# varying parameter :$(parameter)\n#\n")
        # write header line
        hl = "# :$(parameter)"
        for j in 1:length(bands)
            hl = hl * "\tE_$(j)"
        end
        print(f, hl, "\n")
        # write body
        for i in 1:length(values)
            l = "$(values[i])"
            for j in 1:length(bands)
                l = l * "\t$(bands[j][i])"
            end
            print(f, l, "\n")
        end
        # close file
        close(f)
    end
    # plotting
    ylabel("energy")
    xlabel("parameter \"$(parameter)\"")
    xlim(values[1], values[end])
    for b in bands
        plot(values, b, color=color; kwargs...)
    end
    set_parameter!(op, parameter, ba_param, site=site)
end
export show_energy_evolution

function show_energy_evolution_with_parameters(op :: AbstractOperator, parameter::Symbol, values; site=:all, subtract_GS::Bool = false, figsize::Tuple=(8,5))
    figure(figsize=(figsize[1]+6, figsize[2]))

    # annotation
    subplot(122)
    io = IOBuffer()
    show(io, op)
    annotation = String(take!(io)) * "\n\n" * "Parameters:"
    parameter_list = get_parameters(op, site=:all)
    deleteat!(parameter_list, findfirst(p->p==parameter, parameter_list))
    sites = get_sites(basis(op))
    for p in parameter_list
        annotation *= "\n - " * String(p)
        param_values = []
        param_sites  = Vector{Int64}[]
        for s in sites
            pv = get_parameter(op, p, site=s)
            if pv != nothing
                if pv in param_values
                    push!(param_sites[findfirst(e->e==pv, param_values)], s)
                else
                    push!(param_values, pv)
                    push!(param_sites, Int64[s,])
                end
            end
        end
        if length(param_values) == 1
            annotation *= " = "*string(param_values[1])
        else
            annotation *= " = "
            for i in 1:length(param_values)
                annotation *= " (" * string(param_values[i]) * " @sites " * string(param_sites[i][1])
                for s in 2:length(param_sites[i])
                    annotation *= "," * string(param_sites[i][s])
                end
                annotation *= ")"
            end
        end
    end
    gca().text(0.0,1.0,annotation,horizontalalignment="left", verticalalignment="top")
    ylim(-0.05, 1.05)
    setp(gca(), frame_on=false, xticks=(), yticks=())

    # actual plot
    subplot(121)
    bands = [zeros(length(values)) for i in 1:length(basis(op))]
    ba_param = get_parameter(op, parameter, site=site)
    for i in 1:length(values)
        set_parameter!(op, parameter, values[i]; recalculate=true, site=site)
        evals = energies(op)
        for j in 1:length(evals)
            bands[j][i] = evals[j]
        end
        if subtract_GS
            for j in 1:length(evals)
                bands[j][i] -= evals[1]
            end
        end
    end
    # plotting
    ylabel("energy")
    xlabel("parameter \"$(parameter)\"")
    xlim(values[1], values[end])
    for b in bands
        plot(values, b, color="b")
    end
    set_parameter!(op, parameter, ba_param, site=site)

    # tigten the layout
    tight_layout()
end
export show_energy_evolution_with_parameters



function showMPState_evolution(
        op :: AbstractOperator,
        index :: Integer,
        parameter::Symbol,
        values
        ;
        site = :all,
        cutoff :: Real = 0.01,
        digits :: Integer = 3,
        max_contributions :: Integer = 10
    )
    # empty lists for the contributions
    contributions = [zeros(length(values)) for i in 1:length(basis(op))]

    # backup the parameter
    ba_param = get_parameter(op, parameter, site=site)
    # get all contributions
    for i in 1:length(values)
        # set the paramter
        set_parameter!(op, parameter, values[i]; recalculate=true, site=site)
        eigsys = eigensystem(op)
        # obtain the eigenvector
        evec = eigsys[:vectors][index]
        max_weight = sum([abs(v) for v in evec])
        # set the contributions
        for j in 1:length(evec)
            contributions[j][i] = abs(evec[j]) / max_weight
        end
    end
    # reset the parameter to backup
    set_parameter!(op, parameter, ba_param, site=site)

    # calculate the integrated contributions
    contributions_integrated = [sum(contrib)/length(contrib) for contrib in contributions]


    # calculate most important contributions
    contrib_list = [
        (i, contributions[i], contributions_integrated[i]) for i in 1:length(contributions_integrated)
    ]
    contrib_list = sort(contrib_list, by=c->-c[3])

    # obtain basis
    b = basis(op)
    # plot the largest contributions
    ci = 0
    ylabel("contribution [%]")
    xlim(values[1], values[end])
    xlabel("parameter \"$(parameter)\"")
    for c in contrib_list
        # compose the label
        bs = "(+)"
        for sps in b[c[1]].occupation
             bs *= " d_" * summary(b.single_particle_basis[sps])
        end
        sites_bs = unique([b.single_particle_basis[sps].site for sps in b[c[1]].occupation])
        si = length(sites_bs)
        bs *= " |vac>"
        bs_long = bs * "\t ($(round(c[3]*100, digits=2))%)"
        println(bs_long)
        plot(values, c[2].*100, label=bs, color=((si-1)/(s-1), 0, 1-(si-1)/(s-1)))
        ci += 1
        if ci == max_contributions && length(contrib_list) != max_contributions
            println("... (+ $(length(contrib_list)-ci) more)")
            break
        end
    end
    legend()
end
export showMPState_evolution

function showMPState_evolution_grouped(
        op :: AbstractOperator,
        indices :: Vector{<:Integer},
        parameter::Symbol,
        values,
        grouping_fkt :: Function
        ;
        site = :all,
        cutoff :: Real = 0.01,
        digits :: Integer = 3,
        average :: Bool = false
    )

    # make groups
    group_labels = Vector{String}(undef, length(basis(op)))
    group_ids    = zeros(Int64, length(basis(op))) .- 1

    # find out all group labels
    for b in 1:length(group_ids)
        # find out the label and set it
        group_labels[b] = grouping_fkt(basis(op)[b], b, basis(op))
    end
    # produce a unique and sorted label list
    group_labels  = sort(unique(group_labels))
    group_lengths = zeros(Int64, length(group_labels))
    # find out which state has which label
    for b in 1:length(group_ids)
        # find out the label and set it
        group_ids[b] = findfirst(group_labels .== grouping_fkt(basis(op)[b], b, basis(op)))
        group_lengths[group_ids[b]] += 1
    end

    # make a group contribution list
    group_contributions = [zeros(length(values)) for i in 1:length(group_labels)]

    # backup the parameter
    ba_param = get_parameter(op, parameter, site=site)
    # get all contributions
    for i in 1:length(values)
        # set the paramter
        set_parameter!(op, parameter, values[i]; recalculate=true, site=site)
        eigsys = eigensystem(op)
        # obtain the eigenvector
        for index in indices
            evec = eigsys[:vectors][index]
            max_weight = sum([abs(v) for v in evec])
            # get the individual contributions
            contributions = Float64[
                abs(evec[j]) / max_weight for j in 1:length(evec)
            ]
            # bunch to the group contributions
            for j in 1:length(evec)
                group_contributions[group_ids[j]][i] += contributions[j] / length(indices)
            end
        end
    end
    # reset the parameter to backup
    set_parameter!(op, parameter, ba_param, site=site)

    # plot the group contributions
    ci = 0
    if average
        ylabel("group contribution [%] per basis state")
    else
        ylabel("group contribution [%]")
    end
    xlim(values[1], values[end])
    xlabel("parameter \"$(parameter)\"")
    for c in 1:length(group_labels)
        # compose the label
        if average
            plot(values, group_contributions[c].*100 ./ group_lengths[c], label="$(group_labels[c])  ($(group_lengths[c]) total)")
        else
            plot(values, group_contributions[c].*100, label=group_labels[c])
        end
    end
    legend()
end

function showMPState_evolution_grouped(
        op :: AbstractOperator,
        index :: Integer,
        parameter::Symbol,
        values,
        grouping_fkt :: Function
        ;
        site = :all,
        cutoff :: Real = 0.01,
        digits :: Integer = 3,
        average :: Bool = false
    )

    # make groups
    group_labels = Vector{String}(undef, length(basis(op)))
    group_ids    = zeros(Int64, length(basis(op))) .- 1

    # find out all group labels
    for b in 1:length(group_ids)
        # find out the label and set it
        group_labels[b] = grouping_fkt(basis(op)[b], b, basis(op))
    end
    # produce a unique and sorted label list
    group_labels  = sort(unique(group_labels))
    group_lengths = zeros(Int64, length(group_labels))
    # find out which state has which label
    for b in 1:length(group_ids)
        # find out the label and set it
        group_ids[b] = findfirst(group_labels .== grouping_fkt(basis(op)[b], b, basis(op)))
        group_lengths[group_ids[b]] += 1
    end

    # make a group contribution list
    group_contributions = [zeros(length(values)) for i in 1:length(group_labels)]

    # backup the parameter
    ba_param = get_parameter(op, parameter, site=site)
    # get all contributions
    for i in 1:length(values)
        # set the paramter
        set_parameter!(op, parameter, values[i]; recalculate=true, site=site)
        eigsys = eigensystem(op)
        # obtain the eigenvector
        evec = eigsys[:vectors][index]
        max_weight = sum([abs(v) for v in evec])
        # get the individual contributions
        contributions = Float64[
            abs(evec[j]) / max_weight for j in 1:length(evec)
        ]
        # bunch to the group contributions
        for j in 1:length(evec)
            group_contributions[group_ids[j]][i] += contributions[j]
        end
    end
    # reset the parameter to backup
    set_parameter!(op, parameter, ba_param, site=site)

    # plot the group contributions
    ci = 0
    if average
        ylabel("group contribution [%] per basis state")
    else
        ylabel("group contribution [%]")
    end
    xlim(values[1], values[end])
    xlabel("parameter \"$(parameter)\"")
    for c in 1:length(group_labels)
        # compose the label
        if average
            plot(values, group_contributions[c].*100 ./ group_lengths[c], label="$(group_labels[c])  ($(group_lengths[c]) total)")
        else
            plot(values, group_contributions[c].*100, label=group_labels[c])
        end
    end
    legend()
end
export showMPState_evolution_grouped
