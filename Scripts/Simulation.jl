module SM

using Distributions, LinearAlgebra, Interpolations
dir = "/Users/jbusink/Documents/GitHub/Neutronenbron_Simulatie/"
include(dir * "Scripts/Neutron_functions.jl")

Elist = NF.logrange(-5, 7, 1001)
survival_E = NF.cross_section_interpolated(dir * "Data/H1_elastic_scattering_sigma.rtf", Elist) ./
             (NF.cross_section_interpolated(dir * "Data/H1_elastic_scattering_sigma.rtf", Elist) .+ NF.cross_section_interpolated(dir * "Data/H1_fission_sigma.rtf", Elist))
survival_rate = linear_interpolation(Elist, survival_E)

function unit_vector(theta, phi)
    """
    Return a unit vector with polar angles theta and phi.
    """
    a = cos.(theta)
    b = sin.(theta)
    c = cos.(phi)
    d = sin.(phi)
    return Array([b * c, b * d, a])
end

function Maxwell_Boltzmann(T)
    """A Maxwell-Boltzmann distribution of energies is equivalent 
    to a Gamma distributions with k-shape factor 3/2 and scaling of KbT"""
    Kb = 8.617343 * 1e-5 # ev/K
    T = 293
    d = Distributions.Gamma.(3 / 2, Kb .* T)
    E0 = rand(d, 1)[1]
    if E0 < 1e-5
        E0 = 1e-5
        return E0
    end
    return E0
end

function check_survivalHnγ(E)
    if E == 0
        return E
    end
    if E < 1e-5
        E = 1e-5
    end

    if rand() > survival_rate(E)
        E = 0
        return E
    end
    return E
end

function single_collision_random(P_initial, E_initial)
    if E_initial == 0
        # p_lab_new = 0 .* unit_vector(0, 0)
        # Energy_new = 0
        return (0 .* unit_vector(0, 0), 0)
    end

    θ = (rand()) * pi
    ϕ = rand() * 2 * pi

    p_lab = (P_initial ./ norm(P_initial)) * sqrt(2 * E_initial) #if energy is Boltzmann distributed, I correct for it in the momentum.
    p_cm = p_lab / 2.0

    p_cm_new = norm(p_cm) * unit_vector(θ, ϕ)
    p_lab_new = p_cm_new .+ p_cm
    Energy_new = norm(p_lab_new)^2 ./ 2.0

    # println("CHECK", Energy_new - norm(p_lab_new)^2/2 ) # should be zero.
    
    return (p_lab_new, Energy_new)
end

function multiple_collision_random(P_initial, E_initial, n, thermal)
    Energylist, Px, Py, Pz = Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    append!(Energylist, E_initial)
    append!(Px, P_initial[1])
    append!(Py, P_initial[2])
    append!(Pz, P_initial[3])

    i = 0
    ET = 0.03786
    T = 293
    pmagnitude = norm(P_initial)
    E = pmagnitude .^ 2 ./ 2.0
    while i < n
        E = check_survivalHnγ(E)

        
        if thermal == 1
            if E < ET && E != 0
                E = Maxwell_Boltzmann(T)[1]
                # println("check 3 (boltzman)", length(E), " ", E)

            end
        end

        
        P_initial, E = single_collision_random(P_initial, E)

        append!(Energylist, E)
        append!(Px, P_initial[1])
        append!(Py, P_initial[2])
        append!(Pz, P_initial[3])
        i = i .+ 1
    end


    return Energylist, Px, Py, Pz
end

# function save_single_trajectory()
end


# function collision_lab_x(Energy)
#     pmagnitude = sqrt.(2 * Energy)

#     θ = (rand()) * pi
#     ϕ = rand() * 2 * pi
#     p_lab = pmagnitude * Array([1, 0, 0]) #collision in arbitrary x dimension.
#     p_cm = p_lab / 2.0
#     p_cm_new = norm.(p_cm) * unit_vector.(θ, ϕ)
#     p_lab_new = p_cm_new .+ p_cm
#     Energy_new = norm.(p_lab_new)^2 ./ 2.0
#     return (Energy_new, θ, ϕ)
# end

# function collision_lab_z(Energy)
#     pmagnitude = sqrt.(2 * Energy)

#     θ = (rand()) * pi
#     ϕ = rand() * 2 * pi
#     p_lab = pmagnitude * Array([0, 0, 1]) #collision in arbitrary x dimension.
#     p_cm = p_lab / 2.0
#     p_cm_new = norm(p_cm) * unit_vector.(θ, ϕ)
#     p_lab_new = p_cm_new .+ p_cm
#     Energy_new = norm(p_lab_new)^2 ./ 2.0
#     return (Energy_new, θ, ϕ)
# end

# function multiple_collision(E0, n)
#     i = 0
#     ET = 4e-2
#     T = 293
#     Elist = Vector{Float64}()
#     collision_number = Vector{Float64}()
#     push!(Elist, E0)
#     push!(collision_number, 0)
#     Thermalization_number = Vector{Float64}()

#     while i < n
#         E0, θ, ϕ = collision_lab_z(E0)
#         if E0 < ET
#             E0 = Maxwell_Boltzmann(T)[1]
#             # println(i)
#             push!(Thermalization_number, i)
#         end

#         push!(Elist, E0)
#         push!(collision_number, 0)
#         i += 1
#     end
#     return Elist, collision_number, Thermalization_number
# end