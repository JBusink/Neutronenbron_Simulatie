module Neutron_functions

export logrange, MFP, cross_section_H1_SI

function logrange(start, stop, step, base=10)
    return base .^ range(start, stop, step)
end

function MFP(n, sigma)
    """Calculates the MFP, preferably in cm.

    Args:
        n (float): the number density in atoms/cm^3
        sigma (float): the cross section in barn (10^-24 cm^2)/atom

    Returns:
        lambda (float): The mean free path in cm, 
        this is the average distance travelled before a collision (elastic/inelastic)
    """
    return (1 / (n * sigma))
end

function cross_section_H1_SI(E, E0=4.5e4, ET=0.04)
    """   Calculates the cross section of elastic neutron scattering with hydrogen atoms.
    Data from stralingshygiene book.

    Args:
        E (float): Energy of neutron in eV.
        E0 (float, optional): upperbound of the plateau. Defaults to 4.5e4.
        ET (float, optional): lower bound of the plateau, thermal energyscale, default value is kbT. Defaults to 0.04.

    Returns:
        float: cross section value. The value is in units of m^2?
    """
    sigma_plateau = 20.0 #eV
    cross_section = sigma_plateau * (E0^2 / (E0^2 + E^2) + ET^2 / E^2)^0.25
    return cross_section / 1e28
end

function random_distance_trav(Energy)
    λ0 = MFP(n_protons, cross_section_H1_SI(Energy))
    d = Exponential(λ0)
    l = rand(d)
    return l
end

end