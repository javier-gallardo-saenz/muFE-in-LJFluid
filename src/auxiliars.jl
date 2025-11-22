using Parameters
using function_utils


@with_kw struct H_eval_fun_names
    V::String = ""
    dVdr::String = ""
    λV::String = ""
    dλVdr::String = ""
    dλVdλ::String = ""
end


function save_simulation_state(filename::String, p::Vector{<:MVector{d,<:Real}}, q::Vector{<:MVector{d,<:Real}}, m::Real, params::LJ_params, λ::Real, L::SVector{d,R},
    functions::H_eval_fun_names, γ::Real, kbT::Real, δt::Real, step::Int, tag::String)

    open(filename, "w") do io 
        println(io, step)
        println(io, "$(m) $(γ) $(kbT) $(λ) $(δt)")
        println(io, join(L, ""))
        println(io, "$(params.ϵ) $(params.σ)")

        for fun_name in fieldnames(H_eval_fun_names)
            name = getfield(functions, fun_name)
            println(io, "$(fun_name)=$(name)")
        end

        for qi in q
            println(io, join(qi, " "))
        end

        for pi in p
            println(io, join(pi, " "))
        end
    end
end


function load_simulation(filename::String)
    lines = realines(filename)
    idx = 1

    step = parse(Int, lines[idx]); idx += 1
    m, γ, kbT, λ, δt = parse.(Float64, split(lines[idx])); idx += 1
    L = Svector{3, Float64}(parse.(Float64, split(lines[idx]))); idx += 1
    ϵ, σ = parse.(Float64, split(lines[idx])); idx += 1
    params = LJ_params(ϵ, σ)

    functions = H_eval_fun_names()
    for fun_name in fieldnames(H_eval_fun_names)
            _, name = split(line[idx], "=", limit = 2); idx += 1
            setfield!(functions, fun_name, name)
    end

    N = Int((length(lines) - idx + 1)/2)
    q = Vector[]
    p = Vector[]

    for i in 1:N
        push!(q, MVector{3,Float64}(parse.(Float64, split(lines[idx])))); idx += 1
    end

    for i in 1:N
        push!(q, MVector{3,Float64}(parse.(Float64, split(lines[idx])))); idx += 1
    end

    return step, m, γ, kbT, λ, δt, L, params, functions, p, q

end


function resolve_function(name::String)
    return getfield(function_utils, Symbol(name))
end

