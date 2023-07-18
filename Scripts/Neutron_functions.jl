module NF

using CSV, DataFrames, Interpolations, Distributions


#CONSTANTS

Nav = 6.023 * 1e23 #Avagodros number
MH20 = 18.01524 #atomic mass water.
n_protons = (1e6 * 0.99820 / MH20) * 2 * Nav #(= protons/m^3)
n_oxygen = (1e6 * 0.99820 / MH20) * Nav #(= protons/m^3)

#FUNCTIONS

function logrange(start, stop, step, base=10)
    """Logrange converter. Make a logrange list.

    Args:
        start (float): start value of the list.
        stop (float): stop value of the list.
        step (float): number of steps.
        base (float,optional): base of the log range, default is base 10. 
    Returns:
        logrange (list): list of logrange values.
    """
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
        ET (float, optional): lower bound of the plateau, thermal energyscale, default value is kbT. 
        Defaults to 0.04.

    Returns:
        float: cross section value. The value is in units of m^2?
    """
    sigma_plateau = 20.0 #eV
    cross_section = sigma_plateau * (E0^2 / (E0^2 + E^2) + ET^2 / E^2)^0.25
    return cross_section / 1e28
end


function read_cross_section_data(file)
    """Redas cross sectional data file.

    Args:
        file (string): the location and filename of the cross-sectional data.
        
    Returns:
        df (dataframe): return the dataframe with headers "Frame", "Energy_ev" and "cross_section_b"
    """

    df = CSV.read(file, DataFrame, header=["Frame", "Energy_ev", "cross_section_b"], skipto=10, delim="    ",)
    df.cross_section_b = [parse(Float64, strip(df.cross_section_b[i], ['\\', '}'])) for i in 1:1:length(df.cross_section_b)]
    return df
end


function cross_section_interpolated(files, Energy)
    """Interpolates the cross-section data.

    Args:
        files (string): the location of the data file.
        Energy (float): the energy at which the cross section should be evaluated.

    Returns:
        interpolated (float): the cross-section at the provided energy.
    """

    df = read_cross_section_data(files)
    x, y = df.Energy_ev, df.cross_section_b
    interpolated = linear_interpolation(x, y)
    return interpolated(Energy)
end


function random_distance_trav(Energy)
    """generates a random distance given a certain energy. The random distance depends on the mean free path λ_{MFP}, 
    the mean free path is a function of energy. If the energy is zero, we fix the distance to 1e-10 (sub nanometer).
    This avoids division by zero problems.

    Args:
        Energy (float): the (mean of the) random distance is a function of energy, this is the input energy.

    Returns:
        l (float): the random distance.
    """

    if Energy == 0
        return 1e-10
    end
    λ0 = MFP(n_protons, cross_section_H1_SI(Energy))
    d = Exponential(λ0)
    l = rand(d)
    return l
end

function filter_x_y_z(list,cutoff)
    """Filters a generic list to the first element that exceeds the cut-off value.

    Example:
    list = filter_x_y_z([1,2,3,4,5,4,0],5)
    list = [1,2,3,4]

    Args:
        list (list): a list of number.
        cut-off (float): a cut-off value 

    Returns:
        list (list): a new list with that contain the elements up to the cut-off value.
    """
    for i in 1:1:length(list)
        if list[i] > cutoff
            return list[1:i-1]
        end     
    end
    return list
end


end