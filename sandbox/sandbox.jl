#=
GFM_BUSES_sandbox = [31,32,22,29, 36]
GFL_BUSES_sandbox = [30, 34, 32]
GFM_BUSES_origs = Int[]

function plot_apex(csv_file::String)

    df = CSV.read(csv_file, DataFrame)

    # Expected CSV format:
    #
    # list_GFL          bus_number    apex    nadir
    # "30;34;32;33"     30            ...     ...
    #
    # list_GFL is stored as a String in one CSV cell.


    # ---------------------------------------------------------
    # Prepare GFL configuration labels
    # ---------------------------------------------------------

    # Keep a readable label for the x-axis
    #
    # "30;34;32;33"
    # ->
    # "30, 34, 32, 33"

    df.GFL_label = [
        join(strip.(split(config, ";")), ", ")
        for config in df.list_GFL
    ]

    # Count how many buses are in each GFL configuration
    df.GFL_count = [
        length(split(config, ";"))
        for config in df.list_GFL
    ]


    # ---------------------------------------------------------
    # Find unique GFL configurations
    # ---------------------------------------------------------

    groups = unique(
        select(df, [:list_GFL, :GFL_label, :GFL_count])
    )

    # Sort first by number of GFL buses,
    # then by configuration
    sort!(groups, [:GFL_count, :GFL_label])


    # ---------------------------------------------------------
    # Give each GFL configuration its own x-axis position
    # ---------------------------------------------------------

    group_to_x = Dict(
        row.list_GFL => i
        for (i, row) in enumerate(eachrow(groups))
    )

    df.x_position = [
        group_to_x[config]
        for config in df.list_GFL
    ]

    tick_positions = 1:nrow(groups)
    tick_labels = groups.GFL_label


    # ---------------------------------------------------------
    # Create plot
    # ---------------------------------------------------------

    p = plot(
        xlabel = "GFL Bus Configuration",
        ylabel = "Δf [Hz]",

        xticks = (
            tick_positions,
            tick_labels
        ),

        xrotation = 45,

        legend = :outerright,

        bottom_margin = 10Plots.mm,
        top_margin = 10Plots.mm,

        size = (1200, 700)
    )


    # ---------------------------------------------------------
    # Fixed colors based ONLY on generator type
    # ---------------------------------------------------------

    GFL_COLOR = :purple
    SG_COLOR  = :green
    GFM_COLOR = :orange


    # ---------------------------------------------------------
    # Plot data
    # ---------------------------------------------------------

    for bus in sort(unique(df.bus_number))

        bus_df = filter(
            row -> row.bus_number == bus,
            df
        )

        for row in eachrow(bus_df)

            # Convert:
            #
            # "30;34;32;33"
            #
            # into:
            #
            # [30, 34, 32, 33]

            gfl_buses = parse.(
                Int,
                strip.(split(row.list_GFL, ";"))
            )

            
            # Determine generator type for THIS configuration

            if bus in gfl_buses

                point_color = GFL_COLOR

            elseif bus in GFM_BUSES_origs

                point_color = GFM_COLOR

            else

                point_color = SG_COLOR

            end


            scatter!(
                p,
                [row.x_position],
                [row.apex],

                marker = :circle,
                markersize = 6,

                color = point_color,

                label = false
            )

        end
    end


    # ---------------------------------------------------------
    # Add legend entries once
    # ---------------------------------------------------------

    scatter!(
        p,
        [NaN],
        [NaN],

        marker = :circle,
        markersize = 6,

        color = GFL_COLOR,

        label = "GFL"
    )

    scatter!(
        p,
        [NaN],
        [NaN],

        marker = :circle,
        markersize = 6,

        color = SG_COLOR,

        label = "SG"
    )

    scatter!(
        p,
        [NaN],
        [NaN],

        marker = :circle,
        markersize = 6,

        color = GFM_COLOR,

        label = "GFM"
    )


    # ---------------------------------------------------------
    # Add grouping based on number of GFL buses
    # ---------------------------------------------------------

    counts = sort(
        unique(groups.GFL_count)
    )


    # Current y-axis limits
    ylims_current = ylims(p)

    y_range =
        ylims_current[2] -
        ylims_current[1]


    # Leave some room above data for labels
    new_ymax =
        ylims_current[2] +
        0.10 * y_range

    ylims!(
        p,
        ylims_current[1],
        new_ymax
    )


    # Position of "4 GFL", "5 GFL", etc.
    label_y =
        ylims_current[2] +
        0.05 * y_range


    for count in counts

        # Find configurations with this number of GFL buses
        positions = findall(
            groups.GFL_count .== count
        )

        first_x = first(positions)
        last_x  = last(positions)

        midpoint =
            (first_x + last_x) / 2


        # Group title
        annotate!(
            p,
            midpoint,
            label_y,

            text(
                "$count GFL",
                10,
                :center
            )
        )


        # Vertical separator between groups
        if last_x < nrow(groups)

            vline!(
                p,
                [last_x + 0.5],

                color = :gray,
                linestyle = :dash,
                linewidth = 1,

                label = false
            )

        end

    end


    return p

end

function save_max_min_dur_freq_metrics_for_bus_sandbox(test_name::String, list_GFL_busses::Vector{Int}, bus::Int, outdir::String;  fname::Union{Nothing,String} = nothing,)
    # makes a folder in plots called "test_name" if folder doesn't already exist.
    mkpath(outdir)  #plots/test_name
    

    # get the min and max frequency values
    max_val = 2
    min_val = 1

    # format what should be stored
    df_save = DataFrame(list_GFL=join(list_GFL_busses, ";"),  bus_number=bus,  apex=max_val, nadir=min_val)

    # if fname is not yet set, set it to a formatted filename like "incr_GFL_quant_freq_iteration05.csv"
    fname === nothing && (fname = "max_min_freq.csv")
    fpath = joinpath(outdir, fname)

    # check if file exists and has content (>0 bytes)
    file_has_data = isfile(fpath) && filesize(fpath) > 0

    # append without re-writing the header if the file already contains data
    CSV.write(fpath, df_save; append = file_has_data, writeheader = !file_has_data)
    
    return fpath
end


#save_max_min_dur_freq_metrics_for_bus_sandbox("sandbox_test", GFM_BUSES_sandbox, 3, pwd())

plt_sandbox = plot_apex("max_min_freq.csv")
savefig(plt_sandbox, joinpath(pwd(),"_my_sandbox_test_apex_plot.png"))
=#