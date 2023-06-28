import pandas as pd
from pathlib import Path
from scipy.interpolate import interp1d
import numpy as np
from scipy.optimize import minimize

#Interpolates the computed spectra and their uncertainties for given parameters
class Interpolate:
    def __init__(self, mass, channel, final_state, path = None):
        self.mass = mass
        self.channel = channel 
        self.final_state = final_state
        self.path = path

        self._check_input_parameters()

    #Checks if the input parameters are in a proper format
    def _check_input_parameters(self):
        self._check_mass_format()
        self._check_channel_name()
        self._check_final_state_name()

    # Checks if the entry for mass is a number and if it is between the computed bounds.
    def _check_mass_format(self):
        try:
            self.mass = float(self.mass)
            if self.mass < 5.0:
                raise Exception('The dark-matter mass needs to be higher than 5 GeV')
            elif self.mass > 100000:
                raise Exception('The dark-matter mass needs to be lower than 100 TeV')
        except ValueError:
            raise ValueError('The mass entry needs to be a number')

    # Checks if the given final-state name is in the tables.
    def _check_final_state_name(self):
        good_final_state = ['AntiP', 'Ga', 'Nuel', 'Numu', 'Nuta', 'Positrons']
        if self.final_state not in good_final_state:
            raise NameError('{} is not a good final-state entry. Good final-state entries are AntiP, Ga, Nuel, Numu, Nuta, or Positrons'.format(self.final_state))

    # Checks if the given channel name is in the tables.
    def _check_channel_name(self):
        good_channels = ['ee', 'mumu', 'tautau', 'uu', 'dd', 'ss', 'cc', 'bb', 'tt', 'gaga', 'zz', 'ww', 'gg', 'hh']
        if self.channel not in good_channels:
            raise NameError('{} is not a good channel entry. Good channel entries are ee, mumu, tautau, uu, dd, ss, cc, bb, tt, gaga, zz, ww, gg, or hh'.format(self.channel))

    # Fetches the requested data from the appropriate table.
    def _get_data_from_table(self):
        if self.path == None:
            df = pd.read_table('{0}/{1}'.format(Path(__file__).parent.absolute(), 'data/wUncertainty/AtProduction-{}.dat'.format(self.final_state)), sep = '\s\s+', engine='python')
        else:
            df = pd.read_table('{0}/AtProduction-{1}.dat'.format(self.path, self.final_state), sep = '\s\s+', engine='python')
        dm = df['# DM']

        if not '+DHad [{}]'.format(self.channel) in list(df.columns):
            if self.mass in dm.unique():
                x       = list(df[df['# DM'] == self.mass]['x'])
                dNdx    = list(df[df['# DM'] == self.mass]['dN/dx [{}]'.format(self.channel)])
                return x, dNdx

            else:
                mass_lower, mass_upper = self._check_mass_enclosure(dm)
                x_upper       = list(df[df['# DM'] == mass_upper]['x'])
                dNdx_upper    = list(df[df['# DM'] == mass_upper]['dN/dx [{}]'.format(self.channel)])
                x_lower       = list(df[df['# DM'] == mass_lower]['x'])
                dNdx_lower    = list(df[df['# DM'] == mass_lower]['dN/dx [{}]'.format(self.channel)])
                return mass_upper, x_upper, dNdx_upper, mass_lower, x_lower, dNdx_lower

        if self.mass in dm.unique():
            x       = list(df[df['# DM'] == self.mass]['x'])
            dNdx    = list(df[df['# DM'] == self.mass]['dN/dx [{}]'.format(self.channel)])
            pDHad   = list(df[df['# DM'] == self.mass]['+DHad [{}]'.format(self.channel)]) 
            mDHad   = list(df[df['# DM'] == self.mass]['-DHad [{}]'.format(self.channel)]) 
            pDScale = list(df[df['# DM'] == self.mass]['+DScale [{}]'.format(self.channel)]) 
            mDScale = list(df[df['# DM'] == self.mass]['-DScale [{}]'.format(self.channel)]) 
            pDcNS   = list(df[df['# DM'] == self.mass]['+DcNS [{}]'.format(self.channel)]) 
            mDcNS   = list(df[df['# DM'] == self.mass]['-DcNS [{}]'.format(self.channel)])
            return x, dNdx, pDHad, mDHad, pDScale, mDScale, pDcNS, mDcNS

        else:
            mass_lower, mass_upper = self._check_mass_enclosure(dm)
            x_upper       = list(df[df['# DM'] == mass_upper]['x'])
            dNdx_upper    = list(df[df['# DM'] == mass_upper]['dN/dx [{}]'.format(self.channel)])
            pDHad_upper   = list(df[df['# DM'] == mass_upper]['+DHad [{}]'.format(self.channel)]) 
            mDHad_upper   = list(df[df['# DM'] == mass_upper]['-DHad [{}]'.format(self.channel)]) 
            pDScale_upper = list(df[df['# DM'] == mass_upper]['+DScale [{}]'.format(self.channel)]) 
            mDScale_upper = list(df[df['# DM'] == mass_upper]['-DScale [{}]'.format(self.channel)]) 
            pDcNS_upper   = list(df[df['# DM'] == mass_upper]['+DcNS [{}]'.format(self.channel)]) 
            mDcNS_upper   = list(df[df['# DM'] == mass_upper]['-DcNS [{}]'.format(self.channel)])
            x_lower       = list(df[df['# DM'] == mass_lower]['x'])
            dNdx_lower    = list(df[df['# DM'] == mass_lower]['dN/dx [{}]'.format(self.channel)])
            pDHad_lower   = list(df[df['# DM'] == mass_lower]['+DHad [{}]'.format(self.channel)]) 
            mDHad_lower   = list(df[df['# DM'] == mass_lower]['-DHad [{}]'.format(self.channel)]) 
            pDScale_lower = list(df[df['# DM'] == mass_lower]['+DScale [{}]'.format(self.channel)]) 
            mDScale_lower = list(df[df['# DM'] == mass_lower]['-DScale [{}]'.format(self.channel)]) 
            pDcNS_lower   = list(df[df['# DM'] == mass_lower]['+DcNS [{}]'.format(self.channel)]) 
            mDcNS_lower   = list(df[df['# DM'] == mass_lower]['-DcNS [{}]'.format(self.channel)])
            return mass_upper, x_upper, dNdx_upper, pDHad_upper, mDHad_upper, pDScale_upper, mDScale_upper, pDcNS_upper, mDcNS_upper, mass_lower, x_lower, dNdx_lower, pDHad_lower, mDHad_lower, pDScale_lower, mDScale_lower, pDcNS_lower, mDcNS_lower

    # Checks between wich two computed masses the user-specified mass is.
    def _check_mass_enclosure(self, dm):
        masses = dm.unique()
        for m in masses:
            if self.mass < m:
                mass_upper = m
                break
            mass_lower = m
        return mass_lower, mass_upper

    # Computes the linear fit between two points and subsequently returns the corresponding value for a given mass.
    def _linear_function(self, mass_lower, mass_upper, flux_lower, flux_upper):
        a = (flux_lower - flux_upper)/(mass_lower-mass_upper)
        b = (flux_upper*mass_lower - flux_lower*mass_upper)/(mass_lower - mass_upper)
        return a * self.mass + b

    # Call _linear_function for all variations and values of x and returns the computed arrays.
    def _linear_interpolation(self):
        data = self._get_data_from_table()
        if len(data) == 2:
            return data
        elif len(data) == 6:
            dNdx = []
            for i in range(len(data[1])):
                dNdx.append(self._linear_function(data[3], data[0], data[2][i], data[5][i]))
            return data[1], dNdx
        if len(data) == 8:
            return data
        elif len(data) == 18:
            dNdx, pDHad, mDHad, pDScale, mDScale, pDcNS, mDcNS = [], [], [], [], [], [], []
            for i in range(len(data[1])):
                dNdx.append(self._linear_function(data[9], data[0], data[2][i], data[11][i]))
                pDHad.append(self._linear_function(data[9], data[0], data[3][i], data[12][i]))
                mDHad.append(self._linear_function(data[9], data[0], data[4][i], data[13][i]))
                pDScale.append(self._linear_function(data[9], data[0], data[5][i], data[14][i]))
                mDScale.append(self._linear_function(data[9], data[0], data[6][i], data[15][i]))
                pDcNS.append(self._linear_function(data[9], data[0], data[7][i], data[16][i]))
                mDcNS.append(self._linear_function(data[9], data[0], data[8][i], data[17][i]))
            return data[1], dNdx, pDHad, mDHad, pDScale, mDScale, pDcNS, mDcNS

    # A dataframe is being made that will be the final output.
    def _generate_dataframe_output(self, f):
        if len(f) == 2:
            df = pd.DataFrame(data = {'x': f[0], 'dN/dx': f[1]})
            return df
        elif len(f) == 8:
            df = pd.DataFrame(data = {'x': f[0], 'dN/dx': f[1], '+DHad': f[2], '-DHad': f[3], '+DScale': f[4], '-DScale': f[5], '+DcNS': f[6], '-DcNS': f[7]})
            return df

    # The function that will provide the output of the spectrum for a given mass, channel, and final state in a pandas dataframe format.
    def make_spectrum(self):
        f = self._linear_interpolation()
        return self._generate_dataframe_output(f)

#This class diffuses a given spectrum
class Diffuse(Interpolate):
    #Please note that in this class the spectra are given in dN/dE - E by default, as opposed to dN/dx - x
    def __init__(self, mass, channel, final_state, path = None, path_smearing = None):
        self.mass = mass
        self.path = path

        self._check_final_state(final_state)
        self._check_mass()

        self.spectrum_pre = None
        self.spectrum_pre = Interpolate(mass, channel, final_state, path).make_spectrum()
        self._convert_spectrum_pre_to_E()

        self.spectrum_post = None
        self._create_empty_spectrum_post()

        if path_smearing == None:
            self.diffusion_spectrum = pd.read_csv('{0}/{1}'.format(Path(__file__).parent.absolute(), 'Diffusion_values.csv'))
        else:
            self.diffusion_spectrum = pd.read_csv('{}'.format(path_smearing))

    #Checks if the final state is antiprotons. Diffusion only works for these particles. Other final-state particles might follow later.
    def _check_final_state(self, channel):
        if not channel == 'AntiP':
            raise NameError('Diffusion of spectra only works for antiprotons (AntiP)')

    #Checks if the DM mass is not too low, otherwise the results may be unreliable due to a lack of solar modulation
    def _check_mass(self):
        if self.mass < 20:
            raise Warning('Results for the spectrum below 5 GeV may be unreliable')

    #Converts the pre-diffused spectrum from dN/dx -x to dN/dE - E
    def _convert_spectrum_pre_to_E(self):
        self.spectrum_pre['x'] = self.spectrum_pre['x'] * self.mass 
        self.spectrum_pre['dN/dx'] = self.spectrum_pre['dN/dx'] / self.mass
        self.spectrum_pre['+DHad'] = self.spectrum_pre['+DHad'] / self.mass
        self.spectrum_pre['-DHad'] = self.spectrum_pre['-DHad'] / self.mass
        self.spectrum_pre['+DScale'] = self.spectrum_pre['+DScale'] / self.mass
        self.spectrum_pre['-DScale'] = self.spectrum_pre['-DScale'] / self.mass
        self.spectrum_pre['+DcNS'] = self.spectrum_pre['+DcNS'] / self.mass
        self.spectrum_pre['-DcNS'] = self.spectrum_pre['-DcNS'] / self.mass
        self.spectrum_pre = self.spectrum_pre.rename(columns = {'x': 'E', 'dN/dx': 'dN/dE'})

    #Creates an empty dataframe to which the post-diffusion spectrum is to be written
    def _create_empty_spectrum_post(self):
        self.spectrum_post = self.spectrum_pre.copy()
        self.spectrum_post['dN/dE'] = 0.0
        self.spectrum_post['+DHad'] = 0.0
        self.spectrum_post['-DHad'] = 0.0
        self.spectrum_post['+DScale'] = 0.0
        self.spectrum_post['-DScale'] = 0.0
        self.spectrum_post['+DcNS'] = 0.0
        self.spectrum_post['-DcNS'] = 0.0

    #Diffuses the pre-diffusion spectrum and writes the result to the spectrum_post attribute
    def _diffuse_spectrum(self):
        log_spline_dictionary = self._interpolated_diffusion_shapes()
        for bin_index in self.spectrum_post.index:
            E = self.spectrum_post.at[bin_index, 'E']
            if E < 0.0011:
                continue
            else:
                self._diffuse_flux_at_specific_E(E, log_spline_dictionary, bin_index)

    #Computes how the flux at a specific energy of the pre-diffuesed spectrum diffuses
    def _diffuse_flux_at_specific_E(self, E, log_spline_dictionary, bin_index):
        E_low, E_high = self._nearest_computed_diffusion_spectrum(E)
        for index in self.spectrum_post.index:
            bin_diff = self.spectrum_post.at[index, 'E']
            if bin_diff > 0.001:
                diffused_bin_value = self._diffused_bin_value(E, E_low, E_high, bin_diff, log_spline_dictionary) 

                self.spectrum_post.at[index, 'dN/dE'] += diffused_bin_value * self.spectrum_pre.at[bin_index, 'dN/dE']
                self.spectrum_post.at[index, '+DHad'] += diffused_bin_value * self.spectrum_pre.at[bin_index, '+DHad']
                self.spectrum_post.at[index, '-DHad'] += diffused_bin_value * self.spectrum_pre.at[bin_index, '-DHad']
                self.spectrum_post.at[index, '+DScale'] += diffused_bin_value * self.spectrum_pre.at[bin_index, '+DScale']
                self.spectrum_post.at[index, '-DScale'] += diffused_bin_value * self.spectrum_pre.at[bin_index, '-DScale']
                self.spectrum_post.at[index, '+DcNS'] += diffused_bin_value * self.spectrum_pre.at[bin_index, '+DcNS']
                self.spectrum_post.at[index, '-DcNS'] += diffused_bin_value * self.spectrum_pre.at[bin_index, '-DcNS']

    #Given a specific energy E, determines the nearest tabulated diffusion energy
    def _nearest_computed_diffusion_spectrum(self, E):
        values_tabulated_diffusion = self._values_tabulated_diffusion()
        if E in values_tabulated_diffusion:
            return E
        else:
            E_low = max([float(i) for i in values_tabulated_diffusion if float(i) < E])
            E_high = min([float(i) for i in values_tabulated_diffusion if float(i) > E])
            return E_low, E_high

    #Interpolates the value of dN/dE for a specific E at a specific bin edge. Futhermore, computes the integral of the dN/dE - E graph for a given E and start and end points, i.e. the number of diffused particles in that bin.
    def _diffused_bin_value(self, M, M_low, M_high, bin_diff, log_spline_dictionary):
        bin_diff_x = bin_diff / M
        try:
                dNdx_lower_peak = self._convert_log_spline_to_x_base(log_spline_dictionary[str(M_low)], bin_diff_x) 
                dNdx_upper_peak = self._convert_log_spline_to_x_base(log_spline_dictionary[str(M_high)], bin_diff_x) 

                dNdx_peak = dNdx_lower_peak + (M-M_low) / (M_high - M_low) * (dNdx_upper_peak - dNdx_lower_peak)

                return dNdx_peak / M
        except ValueError:
            return 0.0

    #Returns a list of floats of the masses for which the diffusion is computed with DRAGON (the source injection spectrum is a delta peak a the specified energy and the integral in normalised to 1)
    def _values_tabulated_diffusion(self):
        values_tabulated_diffusion = list(self.diffusion_spectrum.columns)
        values_tabulated_diffusion.remove('E')
        return values_tabulated_diffusion

    #Interpolates the diffusion spectra in the dN/dx - x format. Note the deviation from the dN/dE - E default format. Furthermore, due to numerical issues this is done with log[10](y_data). This function should always be used in combination with _convert_log_spline_x_base_to_default
    def _interpolated_diffusion_shapes(self):
        values_tabulated_diffusion = self._values_tabulated_diffusion()
        log_spline_dictionary = {}
        for value in values_tabulated_diffusion:
            M = float(value)
            log_spline_dictionary[value] = interp1d(self.diffusion_spectrum['E'] / M, np.log10(self.diffusion_spectrum[value] * M))
        return log_spline_dictionary

    #Converts the formatting from the interpolated splines that have the y axis in log[10] base and the x and y axis in dN/dx - x to the default standard, i.e. dN/dE - E
    def _convert_log_spline_to_x_base(self, spline, x):
        return 10 ** spline(x)

    #Converts the diffused spectrum from dN/dE - E to dN/dx - x base
    def _convert_spectrum_post_to_x(self):
        self.spectrum_post['E'] = self.spectrum_pre['E'] / self.mass 
        self.spectrum_post['dN/dE'] = self.spectrum_post['dN/dE'] * self.mass
        self.spectrum_post['+DHad'] = self.spectrum_post['+DHad'] * self.mass
        self.spectrum_post['-DHad'] = self.spectrum_post['-DHad'] * self.mass
        self.spectrum_post['+DScale'] = self.spectrum_post['+DScale'] * self.mass
        self.spectrum_post['-DScale'] = self.spectrum_post['-DScale'] * self.mass
        self.spectrum_post['+DcNS'] = self.spectrum_post['+DcNS'] * self.mass
        self.spectrum_post['-DcNS'] = self.spectrum_post['-DcNS'] * self.mass
        self.spectrum_post = self.spectrum_post.rename(columns = {'E': 'x', 'dN/dE': 'dN/dx'})
    
    #Makes and returns the diffused spectrum
    def make_diffused_spectrum(self):
        self._diffuse_spectrum()
        self._convert_spectrum_post_to_x()
        return self.spectrum_post

#This class contains the functions that determine the uncertainties on the mass and <sigma v>
class Uncertainty:
    def __init__(self, mass, channel, final_state, path = None, path_smearing = None, Emin = None, Emax = None):
        self.mass = mass
        self.channel = channel
        self.final_state = final_state
        self.path = path
        self.path_smearing = path_smearing
        self.spectrum = Interpolate(mass, channel, final_state, path).make_spectrum()
        if final_state == 'AntiP':
            self.diffused_spectrum = Diffuse(mass, channel, final_state, path, path_smearing).make_diffused_spectrum()
        self.Emin, self.Emax = self._check_Ebounds(Emin, Emax)

    #Checks if Emin is within the allowed bounds and a float
    def _check_Ebounds(self, Emin, Emax):
        Emin = self._check_Emin(Emin)
        Emax = self._check_Emax(Emax)
        if Emin > Emax:
            raise ValueError('Emax is smaller than Emin')
        return Emin, Emax

    #Check the validity of Emin
    def _check_Emin(self, Emin):
        if not (isinstance(Emin, (float, int))) and (Emin != None):
            raise TypeError('Emin must be a number')
        elif Emin == None:
            Emin = self.spectrum['x'].min()
        else:
            Emin = float(Emin)

        if Emin > self.mass:
            raise ValueError('Emin is larger than the dark matter mass')

        if self.final_state == 'AntiP':
            if Emin == None:
                Emin = 0.002
            elif Emin < 0.002:
                raise ValueError('Emin is outside the allowed bounds')

        return Emin

    #Checks the validity of Emax
    def _check_Emax(self, Emax):
        if not (isinstance(Emax, (float, int))) and (Emax != None):
            raise TypeError('Emax must be a number')
        elif Emax == None:
            Emax = self.mass
        else:
            Emax = float(Emax)

        if Emax > self.mass:
            raise ValueError('Emax is larger than the dark matter mass')

        if (self.final_state == 'AntiP') and (Emax < 0.002):
            raise ValueError('Emin is outside the allowed bounds')

        return Emax

    #Computes the + and - uncertainties on <sigma v>
    def sigmav_uncertainty(self):
        y_scale_p = self._optimizer_plus(self._chi2_height, 1.0)
        y_scale_m = self._optimizer_min(self._chi2_height, 1.0)
        return y_scale_p, y_scale_m
    
    #Computes the + and - uncertainties on the mass
    def mass_uncertainty(self):
        mass_p = self._optimizer_plus(self._chi2_mass, self.mass)
        mass_m = self._optimizer_min(self._chi2_mass, self.mass)
        return mass_p, mass_m

    #Computes the + and - uncertainties on the spectrum
    def spectrum_uncertainty(self):
        spec_p = self._optimizer_plus(self._chi2_spectrum, self.mass)
        spec_m = self._optimizer_min(self._chi2_spectrum, self.mass)
        return spec_p, spec_m

    #Disregards points in the spectrum that have not been diffused
    def _filter_spectrum(self, df, mass):
        dft = df.copy()
        dft = dft[dft['x'] * mass > 0.001]
        dft = dft[dft['dN/dx'] > 0.0]
        return dft

    #Computes the chi^2 value of a scaled version of the spectrum with the original spectrum
    def _chi2_height(self, y_scale):
        if self.final_state == 'AntiP':
            df = self._filter_spectrum(self.diffused_spectrum, self.mass).copy()
        else:
            df = self.spectrum.copy()
        
        df = df[(df['x'] > self.Emin / self.mass) & (df['x'] < self.Emax / self.mass)]

        df['unc+'] = (df['+DHad']**2 + df['+DScale']**2 + df['+DcNS']**2)**0.5
        df['unc-'] = (df['-DHad']**2 + df['-DScale']**2 + df['-DcNS']**2)**0.5

        df['height+'] = y_scale * df['dN/dx'] > df['dN/dx']
        df['height-'] = y_scale * df['dN/dx'] < df['dN/dx']

        df['unc'] = df['unc+'] * df['height+'].astype(int) + df['unc-'] * df['height-'].astype(int)

        df['chi2'] = ((df['dN/dx'] - y_scale * df['dN/dx']) / df['unc']) **2
        return df['chi2'].sum() / len(df.index)

    #Computes chi^2/d.o.f. for a given mass compared to the original spectrum
    def _chi2_spectrum(self, mass):

        if self.final_state == 'AntiP':
            df = self._filter_spectrum(self.diffused_spectrum, self.mass)
            dft = self._filter_spectrum(Diffuse(mass, self.channel, self.final_state, self.path, self.path_smearing).make_diffused_spectrum(), mass)
        else:
            df = self.spectrum.copy()
            dft = Interpolate(mass, self.channel, self.final_state, self.path).make_spectrum()

        df['unc+'] = (df['+DHad']**2 + df['+DScale']**2 + df['+DcNS']**2)**0.5
        df['unc-'] = (df['-DHad']**2 + df['-DScale']**2 + df['-DcNS']**2)**0.5

        f_dndx = interp1d(dft['x'] * mass / self.mass, np.log10(dft['dN/dx'] * self.mass / mass))

        df = df[(df['x'] > dft['x'].min() * mass / self.mass) & (df['x'] < dft['x'].max() * mass / self.mass) & (df['x'] > self.Emin / mass) & (df['x'] < self.Emax / mass)]

        df['dN/dx_t'] = 10**f_dndx(df['x'])

        df['height+'] = df['dN/dx_t'] > df['dN/dx']
        df['height-'] = df['dN/dx_t'] < df['dN/dx']

        df['unc'] = df['unc+'] * df['height+'].astype(int) + df['unc-'] * df['height-'].astype(int)

        df['chi2'] = ((df['dN/dx'] - df['dN/dx_t']) / df['unc']) **2
        return df['chi2'].sum() / len(df.index)

    #Computes the chi^2/d.o.f. for given spectrum compared to the original mass
    def _chi2_mass(self, mass):
        if self.final_state == 'AntiP':
            df = self._filter_spectrum(self.diffused_spectrum, self.mass)
            dft = self._filter_spectrum(Diffuse(mass, self.channel, self.final_state, self.path, self.path_smearing).make_diffused_spectrum(), mass)
        else:
            df = self.spectrum
            dft = Interpolate(mass, self.channel, self.final_state, self.path).make_spectrum()

        dft['unc+'] = (dft['+DHad']**2 + dft['+DScale']**2 + dft['+DcNS']**2)**0.5
        dft['unc-'] = (dft['-DHad']**2 + dft['-DScale']**2 + dft['-DcNS']**2)**0.5

        f_dndx = interp1d(df['x'] * self.mass / mass, np.log10(df['dN/dx'] * mass / self.mass))

        dft = dft[(dft['x'] > df['x'].min() * self.mass / mass) & (dft['x'] < df['x'].max() * self.mass / mass) & (dft['x'] > self.Emin / mass) & (dft['x'] < self.Emax / mass)]

        dft['dN/dx_t'] = 10**f_dndx(dft['x'])

        dft['height+'] = dft['dN/dx_t'] > dft['dN/dx']
        dft['height-'] = dft['dN/dx_t'] < dft['dN/dx']

        dft['unc'] = dft['unc+'] * dft['height+'].astype(int) + dft['unc-'] * dft['height-'].astype(int)

        dft['chi2'] = ((dft['dN/dx'] - dft['dN/dx_t']) / dft['unc']) **2
        return dft['chi2'].sum() / len(dft.index)

    #Optimezes a function for values larger than its initial value
    def _optimizer_plus(self, func, init_val):
        p1 = init_val
        p2 = 1.1 * init_val
        while func(p2) < 1.0:
            p1 = p2
            p2 = 1.1 * p2
        counter = 0
        while abs(p1 - p2) / p1 > 0.0001:
            counter += 1
            if counter > 100:
                raise RuntimeError('Number of iterations has exceeded 100, aborting')
            pt = 0.5 * (p1 + p2)
            if func(pt) == np.inf:
                print('Warning: the chi^2 value is infinity. Most likely there is a region with zero uncertainty. Returning 0 uncertainty') 
                return init_val
            if func(pt) > 1.0:
                p2 = pt
            elif func(pt) < 1.0:
                p1 = pt
            elif func(pt) == 1.0:
                break
        return pt

    #Optimizes a function for values smaller than its initial value
    def _optimizer_min(self, func, init_val):
        p1 = init_val
        p2 = 0.9 * init_val
        while func(p2) < 1.0:
            p1 = p2
            p2 = 0.9 * p2
        counter = 0
        while abs(p1 - p2) / p1 > 0.0001:
            counter += 1
            if counter > 100:
                raise RuntimeError('Number of iterations has exceeded 100, aborting')
            pt = 0.5 * (p1 + p2)
            if func(pt) == np.inf:
                print('Warning: the chi^2 value is infinnity. Most likely there is a region with zero uncertainty. Returning 0 uncertainty')
                return init_val
            if func(pt) > 1.0:
                p2 = pt
            elif func(pt) < 1.0:
                p1 = pt
            elif func(pt) == 1.0:
                break
        return pt


def annihilation_spectrum(mass, channel, final_state, path = None):
    """
    Returns the linearly interpolated data for the given arguments in a pandas dataframe format

    Parameters:     
        mass : float
            The mass of the annihilating particle in GeV.

        channel : str 
            The particles into which the particles annihilate. Accepted arguments are ee, mumu, tautau, uu, dd, ss, cc, bb, tt, gaga, zz, ww, gg, or hh.
        
        final_state : str
            The spectra of final-state particles. Accepted arguments are AntiP, Ga, Nuel, Numu, Nuta, Positrons. These stand for antiprotons, gamma rays, electron neutrinos, muon neutrinos, tau neutrinos, or positrons respectively.

        path : str, optional
            The path to the annihilation spectra tables. The default path is data/wUncertainty. The output contains the QCD-uncertainties by default. In order to solely get the nominal values the path to the tables without the QCD-uncertainties needs to be given.
                        
    Returns:        
        pandas.DataFrame   

    Disclaimer: The interpolated spectra may be unreliable for masses close to resonant channels (e.g. W, Z, etc.) or for TeV-scale masses due to the lack of high-mass data. This code is a work in progress, it is advised that a sanity check is performed for all obtained spectra.
    """
    spectrum = Interpolate(mass, channel, final_state, path)
    return spectrum.make_spectrum()

def diffused_spectrum(mass, channel, final_state, path = None, path_smearing = None):
    """
    Returns the diffused spectrum for the given arguments in a pandas dataframe format

    Parameters:     
        mass : float
            The mass of the annihilating particle in GeV.

        channel : str 
            The particles into which the particles annihilate. Accepted arguments are ee, mumu, tautau, uu, dd, ss, cc, bb, tt, gaga, zz, ww, gg, or hh.
        
        final_state : str
            The spectra of final-state particles. Currently the only accepted argument is AntiP. This stands for antiprotons. Other final-state particles may be added in the future.

        path : str, optional
            The path to the annihilation spectra tables. The default path is data/wUncertainty. The output contains the QCD-uncertainties by default. In order to solely get the nominal values the path to the tables without the QCD-uncertainties needs to be given.

        path_smearing : str, optional
            The path to the interpolated values of the diffusion of various delta peaks by a propagation program.
                        
    Returns:        
        pandas.DataFrame   

    Disclaimer: The propagated spectrum is not properly normalized and should only be used to determine the uncertainties on <sigma v> and the dark matter mass.
    """
    spectrum = Diffuse(mass, channel, final_state, path, path_smearing)
    return spectrum.make_diffused_spectrum()

def mass_uncertainty(mass, channel, final_state, path = None, path_smearing = None, Emin = None, Emax = None):
    """
    Returns the upper and lower bounds respectively of the masses whose spectra lie in the uncertainty bands of the spectrum of the given mass.

    Parameters:     
        mass : float
            The mass of the annihilating particle in GeV.

        channel : str 
            The particles into which the particles annihilate. Accepted arguments are ee, mumu, tautau, uu, dd, ss, cc, bb, tt, gaga, zz, ww, gg, or hh.
        
        final_state : str
            The spectra of final-state particles. Only AntiP gives the uncertainties for a diffused spectrum, all other arguments determine the uncertainties on the annihilation spectrum.

        path : str, optional
            The path to the annihilation spectra tables. The default path is data/wUncertainty. The output contains the QCD-uncertainties by default. In order to solely get the nominal values the path to the tables without the QCD-uncertainties needs to be given.

        path_smearing : str, optional
            The path to the interpolated values of the diffusion of various delta peaks by a propagation program.

        Emin : float, optional
            The minimum energy for which the spectrum is fitted
        
        Emax : float, optional
            The maximum energy for which the spectrum is fitted
                        
    Returns:        
        tuple  

    Disclaimer: The propagated spectrum is not properly normalized and should only be used to determine the uncertainties on <sigma v> and the dark matter mass.
    """
    return Uncertainty(mass, channel, final_state, path, path_smearing, Emin, Emax).mass_uncertainty()

def spectrum_uncertainty(mass, channel, final_state, path = None, path_smearing = None, Emin = None, Emax = None):
    """
    Returns the upper and lower bounds respectively of the masses whose uncertainty enclose the spectrum of the given mass.

    Parameters:     
        mass : float
            The mass of the annihilating particle in GeV.

        channel : str 
            The particles into which the particles annihilate. Accepted arguments are ee, mumu, tautau, uu, dd, ss, cc, bb, tt, gaga, zz, ww, gg, or hh.
        
        final_state : str
            The spectra of final-state particles. Only AntiP gives the uncertainties for a diffused spectrum, all other arguments determine the uncertainties on the annihilation spectrum.

        path : str, optional
            The path to the annihilation spectra tables. The default path is data/wUncertainty. The output contains the QCD-uncertainties by default. In order to solely get the nominal values the path to the tables without the QCD-uncertainties needs to be given.

        path_smearing : str, optional
            The path to the interpolated values of the diffusion of various delta peaks by a propagation program.
        
        Emin : float, optional
            The minimum energy for which the spectrum is fitted
        
        Emax : float, optional
            The maximum energy for which the spectrum is fitted
                        
    Returns:        
        tuple  

    Disclaimer: The propagated spectrum is not properly normalized and should only be used to determine the uncertainties on <sigma v> and the dark matter mass.
    """
    return Uncertainty(mass, channel, final_state, path, path_smearing, Emin, Emax).spectrum_uncertainty()

def sigmav_uncertainty(mass, channel, final_state, path = None, path_smearing = None, Emin = None, Emax = None):
    """
    Returns the upper and lower fractions respectively for which values the spectrum still lies within the uncertainty bounds

    Parameters:     
        mass : float
            The mass of the annihilating particle in GeV.

        channel : str 
            The particles into which the particles annihilate. Accepted arguments are ee, mumu, tautau, uu, dd, ss, cc, bb, tt, gaga, zz, ww, gg, or hh.
        
        final_state : str
            The spectra of final-state particles. Only AntiP gives the uncertainties for a diffused spectrum, all other arguments determine the uncertainties on the annihilation spectrum.

        path : str, optional
            The path to the annihilation spectra tables. The default path is data/wUncertainty. The output contains the QCD-uncertainties by default. In order to solely get the nominal values the path to the tables without the QCD-uncertainties needs to be given.

        path_smearing : str, optional
            The path to the interpolated values of the diffusion of various delta peaks by a propagation program.

        Emin : float, optional
            The minimum energy for which the spectrum is fitted
        
        Emax : float, optional
            The maximum energy for which the spectrum is fitted
                        
    Returns:        
        tuple  

    Disclaimer: The propagated spectrum is not properly normalized and should only be used to determine the uncertainties on <sigma v> and the dark matter mass.
    """
    return Uncertainty(mass, channel, final_state, path, path_smearing, Emin, Emax).sigmav_uncertainty()
