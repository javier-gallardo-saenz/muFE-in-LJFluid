module muFE_LJ

    using StaticArrays
    using Random
    using LinearAlgebra
    using Distances
    using ProgressBars

    include("utils.jl")
    export LJ_params, kB, box_length, T_real, rho_real, mu_real, reduced_time_unit, t_real

    include("function_utils.jl")
    export LJ_pot, LJ_pot_der, λLJ, dλLJ_dr, dλLJ_dλ, λsoft_pot, dλsoft_pot_dr, dλsoft_pot_dλ,
    Hpi, Upi, dHpi_dq!, dHpi_dq_and_dλ!, twoU_and_dHpi_dq!, allU_and_dHpi_dq!, dHpi_dp!

    include("auxiliars.jl")

    include("integrators.jl")
    export BAOAB!, BAOAB_ti!, BAOAB_twoλ!, BAOAB_allλ!, BAOAB_rdf!

    include("fe_methods.jl")
    export widom_method, thermodynamic_integration, collect_fep_bar, collect_mbar, fep, BAR

end