module muFE_LJ

    using StaticArrays
    using Random
    using LinearAlgebra
    using Distances
    using ProgressBars

    include("utils.jl")
    export LJ_params, kB, box_length, T_real, rho_real, mu_real, reduced_time_unit, t_real, autocorr_time

    include("function_utils.jl")
    export LJ_pot, LJ_pot_der, λLJ, dλLJ_dr, dλLJ_dλ, λsoft_pot, dλsoft_pot_dr, dλsoft_pot_dλ, tailcorr_μex_LJ, tailcorr_U_LJ,
    Hpi, Upi, dHpi_dq!, dHpi_dq_and_dλ!, twoU_and_dHpi_dq!, allU_and_dHpi_dq!, dHpi_dp!

    include("auxiliars.jl")

    include("integrators.jl")
    export onestep_BAOAB!, BAOAB_load_and_save!, BAOAB!, BAOAB_ti!, BAOAB_twoλ!, BAOAB_allλ!, BAOAB_rdf!

    include("fe_methods.jl")
    export widom_method, collect_ti, thermodynamic_integration, collect_fep, collect_mbar, fep, BAR, MBAR

end