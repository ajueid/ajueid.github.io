import pandas as pd

# Checks if the entry for mass is a number and if it is between the computed bounds.
def _check_mass_format(mass):
    try:
        mass = float(mass)
        if mass < 5.0:
            raise Exception('The dark-matter mass needs to be higher than 5 GeV')
        elif mass > 100000:
            raise Exception('The dark-matter mass needs to be lower than 100 TeV')
    except ValueError:
        raise ValueError('The mass entry needs to be a number')

# Checks if the given final-state name is in the tables.
def _check_final_state_name(final_state):
    good_final_state = ['AntiP', 'Ga', 'Nuel', 'Numu', 'Nuta', 'Positrons']
    if final_state not in good_final_state:
        raise NameError('{} is not a good final-state entry. Good final-state entries are AntiP, Ga, Nuel, Numu, Nuta, or Positrons'.format(final_state))

# Checks if the given channel name is in the tables.
def _check_channel_name(channel):
    good_channels = ['ee', 'mumu', 'tautau', 'uu', 'dd', 'ss', 'cc', 'bb', 'tt', 'gaga', 'zz', 'ww', 'gg', 'hh']
    if channel not in good_channels:
        raise NameError('{} is not a good channel entry. Good channel entries are ee, mumu, tautau, uu, dd, ss, cc, bb, tt, gaga, zz, ww, gg, or hh')

# Fetches the requested data from the appropriate table.
def _get_data_from_table(mass, channel, final_state, path):
    if path == None:
        df = pd.read_table('data/wUncertainty/AtProduction-{}.dat'.format(final_state), sep = '\s\s+', engine='python')
    else:
        df = pd.read_table('{0}/AtProduction-{1}.dat'.format(path, final_state), sep = '\s\s+', engine='python')
    dm = df['# DM']

    if not '+DHad [{}]'.format(channel) in list(df.columns):
        if mass in dm.unique():
            x       = list(df[df['# DM'] == mass]['x'])
            dNdx    = list(df[df['# DM'] == mass]['dN/dx [{}]'.format(channel)])
            return x, dNdx

        else:
            mass_lower, mass_upper = _check_mass_enclosure(mass, dm)
            x_upper       = list(df[df['# DM'] == mass_upper]['x'])
            dNdx_upper    = list(df[df['# DM'] == mass_upper]['dN/dx [{}]'.format(channel)])
            x_lower       = list(df[df['# DM'] == mass_lower]['x'])
            dNdx_lower    = list(df[df['# DM'] == mass_lower]['dN/dx [{}]'.format(channel)])
            return mass_upper, x_upper, dNdx_upper, mass_lower, x_lower, dNdx_lower

    if mass in dm.unique():
        x       = list(df[df['# DM'] == mass]['x'])
        dNdx    = list(df[df['# DM'] == mass]['dN/dx [{}]'.format(channel)])
        pDHad   = list(df[df['# DM'] == mass]['+DHad [{}]'.format(channel)]) 
        mDHad   = list(df[df['# DM'] == mass]['-DHad [{}]'.format(channel)]) 
        pDScale = list(df[df['# DM'] == mass]['+DScale [{}]'.format(channel)]) 
        mDScale = list(df[df['# DM'] == mass]['-DScale [{}]'.format(channel)]) 
        pDcNS   = list(df[df['# DM'] == mass]['+DcNS [{}]'.format(channel)]) 
        mDcNS   = list(df[df['# DM'] == mass]['-DcNS [{}]'.format(channel)])
        return x, dNdx, pDHad, mDHad, pDScale, mDScale, pDcNS, mDcNS

    else:
        mass_lower, mass_upper = _check_mass_enclosure(mass, dm)
        x_upper       = list(df[df['# DM'] == mass_upper]['x'])
        dNdx_upper    = list(df[df['# DM'] == mass_upper]['dN/dx [{}]'.format(channel)])
        pDHad_upper   = list(df[df['# DM'] == mass_upper]['+DHad [{}]'.format(channel)]) 
        mDHad_upper   = list(df[df['# DM'] == mass_upper]['-DHad [{}]'.format(channel)]) 
        pDScale_upper = list(df[df['# DM'] == mass_upper]['+DScale [{}]'.format(channel)]) 
        mDScale_upper = list(df[df['# DM'] == mass_upper]['-DScale [{}]'.format(channel)]) 
        pDcNS_upper   = list(df[df['# DM'] == mass_upper]['+DcNS [{}]'.format(channel)]) 
        mDcNS_upper   = list(df[df['# DM'] == mass_upper]['-DcNS [{}]'.format(channel)])
        x_lower       = list(df[df['# DM'] == mass_lower]['x'])
        dNdx_lower    = list(df[df['# DM'] == mass_lower]['dN/dx [{}]'.format(channel)])
        pDHad_lower   = list(df[df['# DM'] == mass_lower]['+DHad [{}]'.format(channel)]) 
        mDHad_lower   = list(df[df['# DM'] == mass_lower]['-DHad [{}]'.format(channel)]) 
        pDScale_lower = list(df[df['# DM'] == mass_lower]['+DScale [{}]'.format(channel)]) 
        mDScale_lower = list(df[df['# DM'] == mass_lower]['-DScale [{}]'.format(channel)]) 
        pDcNS_lower   = list(df[df['# DM'] == mass_lower]['+DcNS [{}]'.format(channel)]) 
        mDcNS_lower   = list(df[df['# DM'] == mass_lower]['-DcNS [{}]'.format(channel)])
        return mass_upper, x_upper, dNdx_upper, pDHad_upper, mDHad_upper, pDScale_upper, mDScale_upper, pDcNS_upper, mDcNS_upper, mass_lower, x_lower, dNdx_lower, pDHad_lower, mDHad_lower, pDScale_lower, mDScale_lower, pDcNS_lower, mDcNS_lower

# Checks between wich two computed masses the user-specified mass is.
def _check_mass_enclosure(mass, dm):
    masses = dm.unique()
    for m in masses:
        if mass < m:
            mass_upper = m
            break
        mass_lower = m
    return mass_lower, mass_upper

# Computes the linear fit between two points and subsequently returns the corresponding value for a given mass.
def _linear_function(mass, mass_lower, mass_upper, flux_lower, flux_upper):
    a = (flux_lower - flux_upper)/(mass_lower-mass_upper)
    b = (flux_upper*mass_lower - flux_lower*mass_upper)/(mass_lower - mass_upper)
    return a*mass + b

# Call _linear_function for all variations and values of x and returns the computed arrays.
def _linear_interpolation(mass, channel, final_state, path):
    data = _get_data_from_table(mass, channel, final_state, path)
    if len(data) == 2:
        return data
    elif len(data) == 6:
        dNdx = []
        for i in range(len(data[1])):
            dNdx.append(_linear_function(mass, data[3], data[0], data[2][i], data[5][i]))
        return data[1], dNdx
    if len(data) == 8:
        return data
    elif len(data) == 18:
        dNdx, pDHad, mDHad, pDScale, mDScale, pDcNS, mDcNS = [], [], [], [], [], [], []
        for i in range(len(data[1])):
            dNdx.append(_linear_function(mass, data[9], data[0], data[2][i], data[11][i]))
            pDHad.append(_linear_function(mass, data[9], data[0], data[3][i], data[12][i]))
            mDHad.append(_linear_function(mass, data[9], data[0], data[4][i], data[13][i]))
            pDScale.append(_linear_function(mass, data[9], data[0], data[5][i], data[14][i]))
            mDScale.append(_linear_function(mass, data[9], data[0], data[6][i], data[15][i]))
            pDcNS.append(_linear_function(mass, data[9], data[0], data[7][i], data[16][i]))
            mDcNS.append(_linear_function(mass, data[9], data[0], data[8][i], data[17][i]))
        return data[1], dNdx, pDHad, mDHad, pDScale, mDScale, pDcNS, mDcNS

# A dataframe is being made that will be the final output.
def _generate_dataframe_output(f):
    if len(f) == 2:
        df = pd.DataFrame(data = {'x': f[0], 'dN/dx': f[1]})
        return df
    elif len(f) == 8:
        df = pd.DataFrame(data = {'x': f[0], 'dN/dx': f[1], '+DHad': f[2], '-DHad': f[3], '+DScale': f[4], '-DScale': f[5], '+DcNS': f[6], '-DcNS': f[7]})
        return df

# The function that will provide the output of the spectrum for a given mass, channel, and final state in a pandas dataframe format.
def Data(mass, channel, final_state, path = None):
    """
    Returns the linearly interpolated data for the given arguments in a pandas dataframe format

    Parameters:     mass : float
                        The mass of the annihilating particle in GeV.

                    channel : str 
                        The particles into which the particles annihilate. Accepted arguments are ee, mumu, tautau, uu, dd, ss, cc, bb, tt, gaga, zz, ww, gg, or hh.
                    
                    final_state : str
                        The spectra of final-state particles. Accepted arguments are AntiP, Ga, Nuel, Numu, Nuta, Positrons. These stand for antiprotons, gamma rays, electron neutrinos, muon neutrinos, tau neutrinos, or positrons respectively.

                    path : str, optional
                        The path to the annihilation spectra tables. The default path is qcd-dm.github.io-main/data/wUncertainty. The output contains the QCD-uncertainties by default. In order to solely get the nominal values the path to the tables without the QCD-uncertainties needs to be given.
                        
    Returns:        pandas.DataFrame   

    Disclaimer: The interpolated spectra may be unreliable for masses close to resonant channels (e.g. W, Z, etc.) or for TeV-scale masses due to the lack of high-mass data. This code is a work in progress, it is advised that a sanity check is performed for all obtained spectra.
    """
    _check_mass_format(mass)
    _check_channel_name(channel)
    _check_final_state_name(final_state)
    f = _linear_interpolation(mass, channel, final_state, path)
    return _generate_dataframe_output(f)
