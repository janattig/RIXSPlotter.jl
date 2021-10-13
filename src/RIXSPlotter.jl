module RIXSPlotter

    # using
    using PyPlot
    using RIXSCalculator

    # include files
    include("eigensystem_plotting.jl")
    include("spectrum_plotting.jl")
    include("lab_system_plotting/lab_system_plotting.jl")
    include("lab_system_plotting/momentum_transfer_plotting.jl")
    

end
