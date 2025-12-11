using Parameters
using function_utils

"""
Custom struct to save the names of the functions used for a given MD load_simulation
"""
@with_kw struct H_eval_fun_names
    V::String = ""
    dVdr::String = ""
    λV::String = ""
    dλVdr::String = ""
    dλVdλ::String = ""
end


"""
Save current state of a simulation, together with all necessary info to continue running it in the future
"""
function save_simulation_state(filename::String, p::Vector{MVector{d,Float64}}, q::Vector{MVector{d,Float64}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,R},
    functions::H_eval_fun_names, γ::Real, kbT::Real, δt::Real, step::Int, expectancies::Union{nothing, Float64}=nothing) where {d, R<:Real}

    open(filename, "w") do io 
        println(io, "[STEP]")
        println(io, step)

        println(io, "\n[PARAMS]")
        println(io, "$(m) $(γ) $(kbT) $(λ) $(δt)")

        println(io, "\n[BOX]")
        println(io, join(L, ""))

        println(io, "\n[LJ_PARAMS]")
        println(io, "$(params.ϵ) $(params.σ)")


        println(io, "\n[FUNCTIONS]")
        for fun_name in fieldnames(H_eval_fun_names)
            name = getfield(functions, fun_name)
            println(io, "$(fun_name)=$(name)")
        end

        println(io, "\n[POSITIONS]")
        for qi in q
            println(io, join(qi, " "))
        end

        println(io, "\n[MOMENTUMS]")
        for pi in p
            println(io, join(pi, " "))
        end

        if expectancies !== nothing
            println(io, "\n[EXPECTANCIES]")
            println(io, expectancies)
        end
    end
end


"""
Helper function to read our custom MD simulation files 
"""
function skip_header(header::String, idx::Int)
    s = 
    @assert s == header "Expected header $header but got $s"
end


"""
Read an MD simulation file as given by save_simulation_state
"""
function load_simulation(filename::String)
    lines = realines(filename)
    idx = 1

    header = lines[idx] 
    @assert header == "[STEP]" "Expected header [STEP] but got $header"; idx += 1
    step = parse(Int, lines[idx]); idx += 1

    header = lines[idx] 
    @assert header == "[PARAMS]" "Expected header [PARAMS] but got $header"; idx += 1
    m, γ, kbT, λ, δt = parse.(Float64, split(lines[idx])); idx += 1

    header = lines[idx] 
    @assert header == "[BOX]" "Expected header [BOX] but got $header"; idx += 1
    L = Svector{3, Float64}(parse.(Float64, split(lines[idx]))); idx += 1

    header = lines[idx] 
    @assert header == "[LJ_PARAMS]" "Expected header [LJ_PARAMS] but got $header"; idx += 1
    ϵ, σ = parse.(Float64, split(lines[idx])); idx += 1
    params = LJ_params(ϵ, σ)
    
    header = lines[idx] 
    @assert header == "[FUNCTIONS]" "Expected header [FUNCTIONS] but got $header"; idx += 1
    functions = H_eval_fun_names()
    for fun_name in fieldnames(H_eval_fun_names)
            _, name = split(line[idx], "=", limit = 2); idx += 1
            setfield!(functions, fun_name, name)
    end
    
    header = lines[idx] 
    @assert header == "[POSITIONS]" "Expected header [POSITIONS] but got $header"; idx += 1

    N = Int((length(lines) - idx + 1)/2)
    q = Vector[]
    p = Vector[]

    for i in 1:N
        push!(q, MVector{3,Float64}(parse.(Float64, split(lines[idx])))); idx += 1
    end

    header = lines[idx] 
    @assert header == "[MOMENTUMS]" "Expected header [MOMENTUMS] but got $header"; idx += 1
    for i in 1:N
        push!(q, MVector{3,Float64}(parse.(Float64, split(lines[idx])))); idx += 1
    end

    expectancies = nothing
    if idx ≤ length(lines) && lines[idx] == "[EXPECTANCIES]"
        idx += 1
        expectancies = parse.(Float64, split(lines[idx]))
    end


    return step, m, γ, kbT, λ, δt, L, params, functions, p, q, expectancies

end


"""
Get function with given name from function_utils.jl
"""
function resolve_function(name::String)
    return getfield(function_utils, Symbol(name))
end

