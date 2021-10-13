################################
#  PLOTTING OF THE LAB SYSTEM  #
################################



function plot_lab_system(
            ls :: LabSystem
            ;
            sample_surrounding :: Symbol    = :box,
            global_coordinate_frame :: Bool = true,
            site_surrounding :: Symbol      = :octahedron,
            site_coordinate_frame :: Bool   = true,
            site_index :: Bool              = true,
            site_connections :: Symbol      = :all,
            dpi :: Real                     = 100
        )

    # get random colors for the different site surroundings
    color_list = [(rand(), rand(), rand()) for s in 1:length(ls.sites)]

    # make a new figure
    figure(figsize=(6,6), dpi=dpi)

    # plot all sites
    for s in 1:length(ls.sites)

        # get center coordinate
        cs = ls.sites[s]

        # get coordinates within global frame
        p_x_p = get_in_global_coordinates(ls.sample, cs.position .+ get_in_global_coordinates(cs, [+1,0,0]./2))
        p_x_m = get_in_global_coordinates(ls.sample, cs.position .+ get_in_global_coordinates(cs, [-1,0,0]./2))
        p_y_p = get_in_global_coordinates(ls.sample, cs.position .+ get_in_global_coordinates(cs, [0,+1,0]./2))
        p_y_m = get_in_global_coordinates(ls.sample, cs.position .+ get_in_global_coordinates(cs, [0,-1,0]./2))
        p_z_p = get_in_global_coordinates(ls.sample, cs.position .+ get_in_global_coordinates(cs, [0,0,+1]./2))
        p_z_m = get_in_global_coordinates(ls.sample, cs.position .+ get_in_global_coordinates(cs, [0,0,-1]./2))
        p_c = get_in_global_coordinates(ls.sample, cs.position)

        # scatter a dot at the center
        scatter3D([p_c[1]], [p_c[2]], [p_c[3]], color=color_list[s], alpha=1.0, s=100)

        # maybe plot surrounding of site
        if site_surrounding == :octahedron || site_surrounding == :Octahedron
            # plot octahedron
            plot3D(
                [p_z_p[1], p_x_p[1], p_z_m[1], p_x_m[1], p_z_p[1]],
                [p_z_p[2], p_x_p[2], p_z_m[2], p_x_m[2], p_z_p[2]],
                [p_z_p[3], p_x_p[3], p_z_m[3], p_x_m[3], p_z_p[3]], color=color_list[s], alpha=0.9
            )
            plot3D(
                [p_z_p[1], p_y_p[1], p_z_m[1], p_y_m[1], p_z_p[1]],
                [p_z_p[2], p_y_p[2], p_z_m[2], p_y_m[2], p_z_p[2]],
                [p_z_p[3], p_y_p[3], p_z_m[3], p_y_m[3], p_z_p[3]], color=color_list[s], alpha=0.9
            )
            plot3D(
                [p_y_p[1], p_x_p[1], p_y_m[1], p_x_m[1], p_y_p[1]],
                [p_y_p[2], p_x_p[2], p_y_m[2], p_x_m[2], p_y_p[2]],
                [p_y_p[3], p_x_p[3], p_y_m[3], p_x_m[3], p_y_p[3]], color=color_list[s], alpha=0.9
            )
        elseif site_surrounding == :none || site_surrounding == :None
            # plot no surroundings
        else
            # raise error
            error("site_surrounding=:$(site_surrounding) does not match a recognized keyword")
        end

        # maybe plot the site coordinate frame
        if site_coordinate_frame
            plot3D([p_c[1], p_x_p[1]], [p_c[2], p_x_p[2]], [p_c[3], p_x_p[3]], color="r", lw=1)
            plot3D([p_c[1], p_y_p[1]], [p_c[2], p_y_p[2]], [p_c[3], p_y_p[3]], color="g", lw=1)
            plot3D([p_c[1], p_z_p[1]], [p_c[2], p_z_p[2]], [p_c[3], p_z_p[3]], color="b", lw=1)
            text3D((p_x_p .+ [0.1, 0.0, 0.0])..., "x", color="r", fontsize=10)
            text3D((p_y_p .+ [0.1, 0.0, 0.0])..., "y", color="g", fontsize=10)
            text3D((p_z_p .+ [0.1, 0.0, 0.0])..., "z", color="b", fontsize=10)
        end

        # maybe plot the site index
        if site_index
            # make a text annotation
            text3D((p_c .+ [0.15, 0.0, 0.0])..., "$(s)", color=color_list[s], fontsize=20)
        end

        # show local distortion axis
        if :n in get_parameters(ls.hamiltonian, site=s)
            p_n = get_in_global_coordinates(ls.sample, cs.position .+ get_in_global_coordinates(cs, get_parameter(ls.hamiltonian,:n, site=s)))
            plot3D([p_c[1], p_n[1]], [p_c[2], p_n[2]], [p_c[3], p_n[3]], color="k")
            text3D((p_n .+ [0.1, 0.0, 0.0])..., "n", color="k", fontsize=10)
        end

        # show local magnetic field axis
        if :B_dir in get_parameters(ls.hamiltonian, site=s)
            p_B = get_in_global_coordinates(ls.sample, cs.position .+ get_in_global_coordinates(cs, get_parameter(ls.hamiltonian,:B_dir, site=s)))
            plot3D([p_c[1], p_B[1]], [p_c[2], p_B[2]], [p_c[3], p_B[3]], color="k")
            text3D((p_B .+ [0.1, 0.0, 0.0])..., "B", color="k", fontsize=10)
        end
    end



    # plot connectivity
    if site_connections == :all || site_connections == :All

        for s1 in 1:length(ls.sites)
        for s2 in s1+1:length(ls.sites)

            # get center coordinate
            cs1 = ls.sites[s1]
            cs2 = ls.sites[s2]

            # get coordinates within global frame
            p_c_1 = get_in_global_coordinates(ls.sample, cs1.position)
            p_c_2 = get_in_global_coordinates(ls.sample, cs2.position)

            # plot a connection
            plot3D([p_c_1[1], p_c_2[1]], [p_c_1[2], p_c_2[2]], [p_c_1[3], p_c_2[3]], lw=5, color="k", alpha=0.2)

        end
        end
    end




    # plot sample
    if sample_surrounding == :box
        pos_mean = [
            sum([s.position[i] for s in ls.sites]) for i in 1:3
        ] ./ length(ls.sites)
        box_half_width = [
            maximum([abs(s.position[i]) for s in ls.sites]) + 0.6 for i in 1:3
        ]
        # box goes from mean - halfwidth to mean + halfwidth
        c1 = get_in_global_coordinates(ls.sample, pos_mean .+ (box_half_width.*[-1,-1,-1]))
        c2 = get_in_global_coordinates(ls.sample, pos_mean .+ (box_half_width.*[ 1,-1,-1]))
        c3 = get_in_global_coordinates(ls.sample, pos_mean .+ (box_half_width.*[-1, 1,-1]))
        c4 = get_in_global_coordinates(ls.sample, pos_mean .+ (box_half_width.*[ 1, 1,-1]))
        c5 = get_in_global_coordinates(ls.sample, pos_mean .+ (box_half_width.*[-1,-1, 1]))
        c6 = get_in_global_coordinates(ls.sample, pos_mean .+ (box_half_width.*[ 1,-1, 1]))
        c7 = get_in_global_coordinates(ls.sample, pos_mean .+ (box_half_width.*[-1, 1, 1]))
        c8 = get_in_global_coordinates(ls.sample, pos_mean .+ (box_half_width.*[ 1, 1, 1]))
        faces = [
            [c1,c2,c4,c3,c1],
            [c5,c6,c8,c7,c5],
            [c1,c2,c6,c5,c1],
            [c3,c4,c8,c7,c3],
        ]
        for chain in faces
            plot3D([c[1] for c in chain],[c[2] for c in chain],[c[3] for c in chain], color="k", alpha=0.5, lw=2)
        end
    end

    # set the limits of the axes
    pos_mean = [
        sum([s.position[i] for s in ls.sites]) for i in 1:3
    ] ./ length(ls.sites)
    box_half_width = [
        maximum([abs(s.position[i]) for s in ls.sites]) + 0.6 for i in 1:3
    ]
    sample_extent = max(maximum(box_half_width), minimum(box_half_width))
    xlim(-sample_extent,sample_extent)
    ylim(-sample_extent,sample_extent)
    zlim(-sample_extent,sample_extent)


    # global coordinate system
    if global_coordinate_frame
        plot3D([-1,2] .* sample_extent, [0,0], [0,0], color="r", lw=2)
        plot3D([0,0], [-1,2] .* sample_extent, [0,0], color="g", lw=2)
        plot3D([0,0], [0,0], [-1,2] .* sample_extent, color="b", lw=2)
        text3D([2.15 * sample_extent,0,0]..., "x", color="r", fontsize=15)
        text3D([0,2.15 * sample_extent,0]..., "y", color="g", fontsize=15)
        text3D([0,0,2.15 * sample_extent]..., "z", color="b", fontsize=15)
    end

    # beams
    c_in  = "k"
    c_out = (0.5, 0.5, 0.5)
    plot3D([0,-ls.q_in[1] * sample_extent* 2.0/norm(ls.q_in)], [0,-ls.q_in[2] * sample_extent* 2.0/norm(ls.q_in)], [0,-ls.q_in[3] * sample_extent* 2.0/norm(ls.q_in)], color=c_in,  alpha=1, lw=3)
    text3D((-ls.q_in .* sample_extent* 2.1/norm(ls.q_in) .+ [0,0,0.1])..., "IN", color=c_in, fontsize=15)
    plot3D([0,ls.q_out[1] * sample_extent* 2.0/norm(ls.q_out)], [0,ls.q_out[2] * sample_extent* 2.0/norm(ls.q_out)], [0,ls.q_out[3] * sample_extent* 2.0/norm(ls.q_out)], color=c_out, alpha=1, lw=3)
    text3D((ls.q_out .* sample_extent* 2.1/norm(ls.q_out) .+ [0,0,0.1])..., "OUT", color=c_out, fontsize=15)

    # polarizations
    c_in      = (0.4, 0.0, 0.5)
    c_out     = (0.8, 0.5, 0.9)
    eps_in    = ls.epsilon_in
    eps_out_h = ls.epsilon_out_hor
    eps_out_v = ls.epsilon_out_ver
    # in
    e_in_start = -ls.q_in .* 1.9 * sample_extent/norm(ls.q_in) .- eps_in.*0.2
    e_in_end   = -ls.q_in .* 1.9 * sample_extent/norm(ls.q_in) .+ eps_in.*0.4
    plot3D([e_in_start[1],e_in_end[1]], [e_in_start[2],e_in_end[2]], [e_in_start[3],e_in_end[3]], color=c_in,  alpha=1, lw=3)
    text3D((e_in_end .+ [0,0,0.1])..., "EPS_IN", color=c_in, fontsize=12)
    # out horizontal
    e_out_start = ls.q_out .* 1.9 * sample_extent/norm(ls.q_out) .- eps_out_h.*0.2
    e_out_end   = ls.q_out .* 1.9 * sample_extent/norm(ls.q_out) .+ eps_out_h.*0.4
    plot3D([e_out_start[1],e_out_end[1]], [e_out_start[2],e_out_end[2]], [e_out_start[3],e_out_end[3]], color=c_out,  alpha=1, lw=3)
    text3D((e_out_end .+ [0,0,0.2])..., "EPS_OUT_H", color=c_out, fontsize=12)
    # out vertical
    e_out_start = ls.q_out .* 1.9 * sample_extent/norm(ls.q_out) .- eps_out_v.*0.2
    e_out_end   = ls.q_out .* 1.9 * sample_extent/norm(ls.q_out) .+ eps_out_v.*0.4
    plot3D([e_out_start[1],e_out_end[1]], [e_out_start[2],e_out_end[2]], [e_out_start[3],e_out_end[3]], color=c_out,  alpha=1, lw=3)
    text3D((e_out_end .+ [0,0,0.2])..., "EPS_OUT_V", color=c_out, fontsize=12)

    axis("off")

end