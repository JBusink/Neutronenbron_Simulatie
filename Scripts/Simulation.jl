module SM

using Distributions, LinearAlgebra, Interpolations, Dates
dir = "/Users/jbusink/Documents/GitHub/Neutronenbron_Simulatie/"
include(dir * "Scripts/Neutron_functions.jl")

## BASIC FUNCTIONS

function unit_vector(theta, phi)
    """   Generates an unit vector using a random value of ϕ and θ.

    Args:
        θ (float): a random angle between 0 and π.
        ϕ (float): a random angle between 0 and 2π.

    Returns:
        Array: a random vector (3d) using a spherical coordinate system.
    """

    a = cos.(theta)
    b = sin.(theta)
    c = cos.(phi)
    d = sin.(phi)
    return Array([b * c, b * d, a])
end

function Maxwell_Boltzmann(T)
    """   Generates a random energie using Boltzmann statistics. If the energy is below E<1e-5 eV the energy is fixed at 1e-5 eV. 
    This is a very rare event and is implemented since the interpolated data of the cross-sections is bound to 1e-5 eV.

    Note that I use a Gamma distribution with shape k_{shape} = 3/2 and scaling \theta = KbT to simulate a boltzmann distribution.
    This can be shown to produce the correct distribution, https://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution. 

    Args:
        T (float): the temprature at which we want to draw an energie from the Boltzmann distribution.

    Returns:
        E (float): the new energie.
    """
    Kb = 8.617343 * 1e-5 # ev/K
    T = 293
    d = Distributions.Gamma.(3 / 2, Kb .* T)
    E = rand(d, 1)[1]
    if E < 1e-5
        E = 1e-5
        return E
    end
    return E
end


# Simple functions that interpolate the data and interpolates the survival rate as a function of E. 
# The survival_rate function needs as Args: E (float) and return the probability that it survives.
Elist = NF.logrange(-5, 7, 1001)
survival_E = NF.cross_section_interpolated(dir * "Data/H1_elastic_scattering_sigma.rtf", Elist) ./
             (NF.cross_section_interpolated(dir * "Data/H1_elastic_scattering_sigma.rtf", Elist) .+ NF.cross_section_interpolated(dir * "Data/H1_fission_sigma.rtf", Elist))
survival_rate = linear_interpolation(Elist, survival_E)


function check_survivalHnγ(E)
    """  Checks if the particle survives. If the energy of the particle is zero, the functions returns E = 0 (nothing happens).
    If the  energy is below 1e-5 eV, the energy is fixed at 1e-5 eV (to avoid interpolation errors). 
    Finally, we draw a random number x = rand(). If the random number is larger than the survival rate x > survival_rate(E)
    the neutron undergoes a fission reaction (i.e. an elastic collision), hence the energy is zero and the neutron is removed.

    If x < survival_rate(E) the neutron is not removed removed and the energy is unchanged.

    Args:
        E (float): the energy (eV) of the neutron.

    Returns:
        E (float): the energy (eV) of the neutron.
    """
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

#Here I calculate the fluxes. When the particle at time t_i+1 is at a position x>r and t_i x < r the object we add a (right-sided) flux + 1. \
#The same is true for the left-sided flux, but the signs are swapped. I had to add position 0 to the list.
function flux_single_run_right(list1, cutoff)
    """Calculates the flux (from in to outward). Every time a neutrons passes an imaginary barrier at x = cut-off 
    we add 1 to a counter. A neutron can pass multiple times the same barrier. In this case 
    we calculate only the flux from inward to outward.

    Args:
        list1 (list): a list of radial positions.
        cut-off (float): the location of the barrier (in cm).

    Returns:
        flux (int): the flux at a given barrier.
    """
    list = [0; list1]
    flux = 0
    for i in 1:1:(length(list)-1)
        j = i + 1
        if list[j] .> cutoff && list[i] .< cutoff
            flux += 1
        end
    end
    return flux
end

function flux_single_run_left(list1, cutoff)
    """Calculates the flux (from out to inward). Every time a neutrons passes an imaginary barrier at x = cut-off 
    we add 1 to a counter. A neutron can pass multiple times the same barrier. In this case 
    we calculate only the flux from outward to inward.

    Args:
        list1 (list): a list of radial positions.
        cut-off (float): the location of the barrier (in cm).

    Returns:
        flux (int): the flux at a given barrier.
    """
    list = [0; list1]
    flux = 0
    for i in 1:1:(length(list)-1)
        j = i + 1
        if list[j] .< cutoff && list[i] .> cutoff
            flux += 1
        end
    end
    return flux
end

function flux_total(list, cutoff)
    """Calculates the total flux (left + right)

    Args:
        list (list): a list of radial positions.
        cut-off (float): the location of the barrier (in cm).

    Returns:
        flux (int): the total flux at a given barrier.
    """
    # """Returns the total flux  flux"""
    flux_left = flux_single_run_left(list, cutoff)
    flux_right = flux_single_run_right(list, cutoff)
    return flux_left + flux_right
end

function flux_netto(list, cutoff)
    """Calculates the total flux (right - left)

    Args:
        list (list): a list of radial positions.
        cut-off (float): the location of the barrier (in cm).

    Returns:
        flux (int): the netto flux at a given barrier.
    """
    # """Returns the netto flux"""
    flux_left = flux_single_run_left(list, cutoff)
    flux_right = flux_single_run_right(list, cutoff)
    return flux_right - flux_left
end



# COLLISION FUNCTIONS

function single_collision_random(P_initial, E_initial)
    """  Simulates a single collision starting using a random starting angle θ and ϕ. The initial momentum vector is random and normalized.
    The momentum vector is translated into a CM momentum vector where the direction is randomized. The new CM momentum vector is transformed
    back to the lab frame by adding halve the original momemtum vector. 


    Args:
        P_initial (vector): a random 3d vector.
        E_initial (float): the initial energy of the neutron. 

    Returns:
        p_lab_new (vector): the new random momentum vector.
        Energy_new (float): the new energy of the momentum vector.
    """

    if E_initial == 0
        # p_lab_new = 0 .* unit_vector(0, 0)
        # Energy_new = 0
        return (0 .* unit_vector(0, 0), 0)
    end

    θ = acos(2 * rand() - 1)
    ϕ = rand() * 2 * pi

    p_lab = (P_initial ./ norm(P_initial)) * sqrt(2 * E_initial) #if energy is Boltzmann distributed, I correct for it in the momentum.
    p_cm = p_lab / 2.0

    p_cm_new = norm(p_cm) * unit_vector(θ, ϕ)
    p_lab_new = p_cm_new .+ p_cm
    Energy_new = norm(p_lab_new)^2 ./ 2.0

    return (p_lab_new, Energy_new)
end


function multiple_collision_random(P_initial, E_initial, n, thermal)
    """ Simulates multiple random collisions of a neutron. The simulation takes into account a thermal cut-off, below the thermal cut-off
    value the energy is randomized using a Boltzmann distribution with temperature T. Furthermore, the particle can undergo a fission reaction,
    when this happens, the energy and momentum of the neutron is set at 0.

    Args:
        P_initial (vector): the initial momentum vector of the neutron.
        E_initial (float): the initial energy of the neutron.
        n (Integer): number of collisions that we simulate.
        thermal (Integer): if thermal == 1 we use the thermal cut-off (Boltzmann energy distribution at E<ET), otherwise we dont.
    Returns:
        Energylist (list): the energy (eV) of the neutron at each collision.
        Px (list): the x momentum of the neutron at each collisions.
        Py (list): the y momentum of the neutron at each collisions.
        Pz (list): the z momentum of the neutron at each collisions.
    """

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


function Simulate_multiple_collisions(n_particles, n_collisions)
    """ Simulates multiple random collisions of a n neutrons. The simulation takes into account a thermal cut-off, below the thermal cut-off
    value the energy is randomized using a Boltzmann distribution with temperature T. Furthermore, the particle can undergo a fission reaction,
    when this happens, the energy and momentum of the neutron is set at 0.

    The SM.multiple_collision_random() function return a vector P and the energy E during each collision. The position of the neutron
    is calculated by drawing a random distance d in cm and calculating the differential Δx ,Δy and Δz. The sign and value of the differential
    depends on the relative size of each vector, hence we divide by the norm. Finally, we calculate the full trajectory by cummulative summation
    of the differentials.


    Args:
        n_particles (int): the number of neutrons to simulate.
        n_collisions (int): the number of collisions of each neutron.   

    Returns:
        Xtotal (list): the x-position of the neutron at each collision. The list contains n-sublists for each neutron.
        Ytotal (list): the y-position of the neutron at each collision. The list contains n-sublists for each neutron.
        Ztotal (list): the z-position of the neutron at each collision. The list contains n-sublists for each neutron.
        Etotal (list): the energy (eV) of the neutron at each collision. The list contains n-sublists for each neutron.
        Pxtotal (list): the x momentum of the neutron at each collisions. The list contains n-sublists for each neutron.
        Pytotal (list): the y momentum of the neutron at each collisions. The list contains n-sublists for each neutron.
        Pztotal (list): the z momentum of the neutron at each collisions. The list contains n-sublists for each neutron.
    """

    Etotal = []
    Pxtotal, Pytotal, Pztotal = [], [], []
    Xtotal, Ytotal, Ztotal = [], [], []
    for i in 1:1:n_particles
        if i % 10000 == 0
            println(i)
            println("Free momory: ", Sys.free_memory() / 2^30)
            println("Local Time: ", Dates.format(now(), "HH:MM:SS"))
        end
        E0, θ, ϕ = (rand() * 4 .+ 3) .* 1e6, acos(2 * rand() - 1), rand() * 2 * π
        P_initial = sqrt.(2 .* E0) .* SM.unit_vector(θ, ϕ)
        E, Px, Py, Pz = SM.multiple_collision_random(P_initial, E0, n_collisions, 1)

        normlist = [norm(permutedims(Array([Px Py Pz]))[1:3, i]) for i in 1:1:length(Pz)] .+ 1e-20
        Δx = Px .* NF.random_distance_trav.(E) .* 100 ./ normlist
        Δy = Py .* NF.random_distance_trav.(E) .* 100 ./ normlist
        Δz = Pz .* NF.random_distance_trav.(E) .* 100 ./ normlist
        x, y, z = cumsum(Δx), cumsum(Δy), cumsum(Δz)

        push!(Etotal, E)
        push!(Pxtotal, Px)
        push!(Pytotal, Py)
        push!(Pztotal, Pz)
        push!(Xtotal, x)
        push!(Ytotal, y)
        push!(Ztotal, z)
    end
    return (Xtotal, Ytotal, Ztotal, Etotal, Pxtotal, Pytotal, Pztotal)
end

end
