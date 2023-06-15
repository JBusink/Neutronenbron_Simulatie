module SM

using Distributions, LinearAlgebra

function unit_vector(theta, phi)
    """
    Return a unit vector with polar angles theta and phi.
    """
    a = cos(theta)
    b = sin(theta)
    c = cos(phi)
    d = sin(phi)
    return Array([b * c, b * d, a])
end

function collision_lab_x(Energy)
    pmagnitude = sqrt.(2 * Energy)

    θ = (rand()) * pi
    ϕ = rand() * 2 * pi
    p_lab = pmagnitude * Array([1, 0, 0]) #collision in arbitrary x dimension.
    p_cm = p_lab / 2.0
    p_cm_new = norm(p_cm) * unit_vector(θ, ϕ)
    p_lab_new = p_cm_new .+ p_cm
    Energy_new = norm(p_lab_new)^2 ./ 2.0
    return (Energy_new, θ, ϕ)
end

function collision_lab_z(Energy)
    pmagnitude = sqrt.(2 * Energy)

    θ = (rand()) * pi
    ϕ = rand() * 2 * pi
    p_lab = pmagnitude * Array([0, 0, 1]) #collision in arbitrary x dimension.
    p_cm = p_lab / 2.0
    p_cm_new = norm(p_cm) * unit_vector(θ, ϕ)
    p_lab_new = p_cm_new .+ p_cm
    Energy_new = norm(p_lab_new)^2 ./ 2.0
    return (Energy_new, θ, ϕ)
end

function single_collision_random(P_initial, E_initial)
    pmagnitude = norm(P_initial)
    E = pmagnitude .^ 2 ./ 2.0

    θ = (rand()) * pi
    ϕ = rand() * 2 * pi
    p_lab = P_initial
    p_cm = p_lab / 2.0

    p_cm_new = norm(p_cm) * unit_vector(θ, ϕ)
    p_lab_new = p_cm_new .+ p_cm
    Energy_new = norm(p_lab_new)^2 ./ 2.0
    return (p_lab_new, Energy_new)
end

function Maxwell_Boltzmann(T)
    """A Maxwell-Boltzmann distribution of energies is equivalent 
    to a Gamma distributions with k-shape factor 3/2 and scaling of KbT"""
    Kb = 8.617343 * 1e-5 # ev/K
    T = 293
    d = Distributions.Gamma(3 / 2, Kb .* T)
    E0 = rand(d, 1)
    return E0
end

function multiple_collision(E0, n)
    i = 0
    ET = 4e-2
    T = 293
    Elist = Vector{Float64}()
    collision_number = Vector{Float64}()
    push!(Elist, E0)
    push!(collision_number, 0)
    Thermalization_number = Vector{Float64}()

    while i < n
        E0, θ, ϕ = collision_lab_z(E0)
        if E0 < ET
            E0 = Maxwell_Boltzmann(T)[1]
            # println(i)
            push!(Thermalization_number, i)
        end

        push!(Elist, E0)
        push!(collision_number, 0)
        i += 1
    end
    return Elist, collision_number, Thermalization_number
end

function multiple_collision_random(P_initial, E_initial, n, thermal)
    Energylist, Px, Py, Pz = Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    collision_number = Vector{Float64}()
    # println(FPT)
    append!(collision_number, 0)
    append!(Energylist, E_initial)
    append!(Px, P_initial[1])
    append!(Py, P_initial[2])
    append!(Pz, P_initial[3])

    i = 0
    ET = 4e-2
    T = 293
    pmagnitude = norm(P_initial)
    E = pmagnitude .^ 2 ./ 2.0

    Thermal_list1 = Vector{Float64}()
    Thermal_list2 = Vector{Float64}()
    Thermal_list3 = Vector{Float64}()
    Thermal_list4 = Vector{Float64}()

    while i < n
        P_initial, E = single_collision_random(P_initial, E)

        if thermal == 1
            if E < ET
                E = Maxwell_Boltzmann(T)[1]
            end
        end

        append!(Energylist, E)
        append!(Px, P_initial[1])
        append!(Py, P_initial[2])
        append!(Pz, P_initial[3])
        append!(collision_number, i)
        i = i .+ 1
    end


    return Energylist, Px, Py, Pz
end

# function save_single_trajectory()
end