from .QCD_Spectra import Interpolate
from sys import argv
from numpy import cosh, sinh, arccosh, log, floor, log10, sqrt, sign
from numpy.random import rand
import matplotlib.pyplot as plt 
from random import uniform
from bisect import bisect
from tqdm import tqdm
import pandas as pd
from matplotlib.style import use
from matplotlib import rc

use(['science', 'grid', 'notebook'])
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#Gets the input for the spectrum
class Input():
    #Makes the input, can either be done directly by the user, of via an input file.
    def __init__(self):
        try: 
            with open(argv[1]) as fopen:
                f = fopen.read().splitlines()
                self._check_file_input(f)
                self.ann = bool(1 - int(f[0]))
                self.final_state = f[1]
                self.XX = bool(1 - int(f[2]))
                self.mdm = float(f[3])
                if not self.ann:
                    self.mdm = self.mdm * 0.5 
                self.mX = float(f[4].split()[0])
                if len(f[4].split()) == 2:
                    self.pol = str(f[4].split()[1])
                else:
                    self.pol = None
                self._check_ann_kinematics()
                self.decay_mode = self._file_decay_modes(f[5])
                self.nsamples = int(f[6])
                self.csv_path = f[7]
        except (IndexError, FileNotFoundError):
            self.ann = self._get_ann_input()
            self.final_state = self._get_final_state()
            self.XX = self._get_input_XX()
            self.mdm, self.mX, self.pol = self._get_input_masses()
            if not self.ann:
                self.mdm = self.mdm * 0.5
            self._check_ann_kinematics()
            self.decay_mode = self._get_decay_modes()
            self.nsamples, self.csv_path = self._get_input_samples()

    #Checks if the input file has the correct format
    def _check_file_input(self, f):
        if f[1] in ['Nuel', 'Numu', 'Nuta']:
            fs = 'v'
        else:
            fs = f[1].lower()
        try: 
            int(f[0])
        except ValueError:
            raise ValueError('Specification of DM decay or annihilation is wrong, must be 1 for decay or 2 for annihilation.')
        try:
            int(f[2])
        except ValueError:
            if self.ann:
                raise ValueError(f'Specification of process is wrong, must be either 1 for DM DM > {fs} X or 2 for DM DM > XX. Aborting.')
            else:
                raise ValueError(f'Specification of process is wrong, must be either 1 for DM > {fs} X or 2 for DM > XX. Aborting.')
        if not int(f[2]) in [1,2]:
            raise ValueError('Specification of process is wrong. Aborting')
        if not f[1] in ['Nuel', 'Numu', 'Nuta', 'Ga']:
            raise ValueError(f'Improper final state. Good final states are Nuel, Numu, Nuta, or Ga not {f[1]}. Aborting')
        try:
            float(f[3])
            mX = f[4].split()[0]
            float(mX)
        except ValueError:
            raise ValueError('Either the DM mass or X mass is not a float. Aborting')
        if (float(f[3]) < 0 ) or (float(f[4].split()[0]) < 0):
            raise ValueError('The DM mass or X mass is smaller than 0. Aborting')
        f4 = f[4].split()
        if len(f4) == 2:
            if not f4[1] in ['-1', '-0.5' ,'0', '0.5', '1', 'T', 'L']:
                raise ValueError('Format of polarization of X is wrong. Must be -1, -0.5, 0, 0.5, 1, L or T. Aborting')
        try:
            int(f[6])
        except ValueError:
            raise ValueError('Number of samples is not an integer')

    #Gets the user input if the process is DM decay or DM annihilation.
    def _get_ann_input(self):
        ann_input = input('Dark Matter decay (1) or annihilation (2):\n')
        while not int(ann_input) in [1,2]:
            ann_input = input('input must be either 1 for DM decay, or 2 for DM annihilation:\n')
        return bool(1-int(ann_input))

    #Gets the input if the DM annihilation channel is X X or X final state.
    def _get_input_XX(self):
        if self.final_state in ['Nuel', 'Numu', 'Nuta']:
            fs = 'v'
        else:
            fs = self.final_state.lower()
        if self.ann:
            XX_input = input(f'Dark Matter annihilation channel:\n1: DM DM > {fs} X\n2: DM DM > X X\n')
        else:
            XX_input = input(f'Dark Matter decay channel:\n1: DM > {fs} X\n2: DM > X X\n')
        check = self._check_XX_input(XX_input)
        while not check:
            XX_input = input(f'Good answers are 1 or 2, not {XX_input}. Please provide new input:\n')
            check = self._check_XX_input(XX_input)
        return bool(1 - int(XX_input))

    #Checks if the input for the decay mode is 1 or 2
    def _check_XX_input(self, inp):
        try:
            tmp = int(inp)
        except ValueError:
            return False
        if tmp == 1 or tmp == 2:
            return True
        else:
            return False

    #Gets the final state
    def _get_final_state(self):
        final_state = input('Final state particle:\nNuel: Electron Neutrino\nNumu: Muon neutrino\nNuta: Tau Neutrino\nGa:   Photon\n')

        while not final_state in ['Nuel', 'Numu', 'Nuta', 'Ga']:
            final_state = input(f'Final state must be Nuel, Numu, Nuta, or Ga, not {final_state}.\n')
        return final_state

    #Gets the branching ratios and decay modes of the BSM particle X.
    def _get_decay_modes(self):
        br = 0
        X_decay_modes = []
        decay_mode = input('Decay modes of particle X (e.g. 1.0,v,h, or 0.5,v,Z,0.33):\n')
        check = self._check_decay_mode(decay_mode)
        while not check:
            decay_mode = input('Decay mode not in a proper format, a good format for example is 0.5,h,z.\n')
            check = self._check_decay_mode(decay_mode)

        comp = decay_mode.split(',')
        if len(comp) == 3:
            X_decay_modes.append(Decay_modes(float(comp[0]), str(comp[1]).lower(), str(comp[2]).lower(), self.mX, self.final_state))
        elif len(comp) == 4:
            X_decay_modes.append(Decay_modes(float(comp[0]), str(comp[1]).lower(), str(comp[2]).lower(), self.mX, self.final_state, float(comp[3])))
        br += float(comp[0])
        while br < 1.0:
            decay_mode = input('Total branching ratio is less than one, please provide additional decay modes:\n')
            check = self._check_decay_mode(decay_mode)
            while not check:
                decay_mode = input('Decay mode not in a proper format, a good format for example is 0.5 h z.\n')
                check = self._check_decay_mode(decay_mode)

            comp = decay_mode.split(',')
            if len(comp) == 3:
                X_decay_modes.append(Decay_modes(float(comp[0]), str(comp[1]).lower(), str(comp[2]).lower(), self.mX, self.final_state))
            elif len(comp) == 4:
                X_decay_modes.append(Decay_modes(float(comp[0]), str(comp[1]).lower(), str(comp[2]).lower(), self.mX, self.final_state, float(comp[3])))
                br += float(decay_mode.split()[0])
        if br > 1.0:
            raise ValueError('Branching ratio is larger than 1.0: Aborting.')
        return X_decay_modes

    #Checks if the decay mode is in a proper format
    def _check_decay_mode(self, inp):
        SM_particles = ['u', 'd', 's', 'c', 'b', 't', 'e','mu', 'tau', 'w', 'z', 'h', 'g', 'ga','v']
        if len(inp.split(',')) != 3 and len(inp.split(',')) != 4:
            return False
        try:
            float(inp.split(',')[0])
        except ValueError:
            return False
        if (float(inp.split(',')[0]) < 0):
            return False
        elif not (inp.split(',')[1] in SM_particles) or not (inp.split(',')[2] in SM_particles):
            return False
        if len(inp.split(',')) == 4:
            try:
                float(inp.split(',')[3])
            except ValueError:
                return False
            if abs(float(inp.split(','))) > 1.0:
                return False
        else:
            return True

    #Gets the decay modes from the input file    
    def _file_decay_modes(self, modes):
        entries = modes.split() 
        comps = [entries[i].split(',') for i in range(len(entries))]
        for comp in comps:
            if len(comp) != 3 and len(comp) != 4:
                raise ValueError('Branching modes of X not properly defined. Aborting')
        X_decay_modes = []
        if sum(float(comps[i][0]) for i in range(len(comps))) != 1.0:
            raise ValueError('Branching ratio of X does not total 1.0. Aborting')
        for comp in comps:
            if len(comp) == 3:
                X_decay_modes.append(Decay_modes(float(comp[0]), str(comp[1]).lower(), str(comp[2]).lower(), self.mX, self.final_state))
            elif len(comp) == 4:
                X_decay_modes.append(Decay_modes(float(comp[0]), str(comp[1]).lower(), str(comp[2]).lower(), self.mX, self.final_state, float(comp[3])))
        return X_decay_modes

    #Gets the masses of the DM particle and particle X, all SM masses are tabulated.
    def _get_input_masses(self):
        mdm = input('Please provide the mass of the Dark Matter particle in GeV:\n')
        check = self._check_masses(mdm)
        while not check:
            mdm = input('Mass is either not a number or is negative, please provide a new Dark Matter mass:\n')
            check = self._check_masses(mdm)
        mX = input('Please provide the mass of particle X in GeV:\n')
        check = self._check_masses(mX)
        while not check:
            mX = input('Mass is either not a number or is negative, please provide a new particle X mass:\n')
            check = self._check_masses(mX)
        pol = input('Please provide the polarization of X (optional):\n')
        if pol == '\n':
            pol = None
        while not pol in ['-1', '-0.5', '0', '0.5' '1', 'L', 'T' , '\n']:
            pol = input('Polarization input must be -1, -0.5, 0, 0.5, 1, L, T, or an enter. Please provide a new polarization:\n')
        return float(mdm), float(mX), pol

    #Check if the input given for mass is a float and not negative
    def _check_masses(self, inp):
        try:
            tmp = float(inp)
        except ValueError:
            return False
        if tmp < 0:
            return False
        else:
            return True

    #Checks if the annihilation channel and particle masses are allowed.
    def _check_ann_kinematics(self):
        if self.XX:
            if self.mX > self.mdm:
                if self.ann:
                    raise Exception('The mass of particle X is larger than the Dark Matter mass; two Dark Matter particles annihilating into an XX pair is thus kinematically forbidden. Aborting')
                else:
                    raise Exception('The mass of particle X is larger than half the Dark Matter mass; a Dark Matter particle decaying into an XX pair is thus kinematically forbidden. Aborting')
        else:
            if self.mX > 2 * self.mdm:
                if self.ann:
                    raise Exception('The mass of particle X is more than twice the Dark Matter mass; two Dark Matter particles annihilating into a final state and an X is thus kinematically forbidden. Aborting.')
                else:
                    raise Exception('The mass of particle X is more than the Dark Matter mass; a Dark Matter particle decaying into a final state and X is thus kinematically forbidden. Aborting.')

    #Gets the number of required points from the spectrum and the path to the csv file where it is to be saved.
    def _get_input_samples(self):
        nsamples = int(input('Please provide the number of samples from the spectrum:\n'))
        path = str(input("Path to the csv file to which the samples should be saved (or enter 'plot' to plot the samples):\n"))
        return nsamples, path

#Samples the spectrum
class Sample():
    #Gets the parameters from the Input class.
    def __init__(self):
        inp = Input()
        self.ann = inp.ann
        self.XX = inp.XX
        self.final_state = inp.final_state
        self.mdm = inp.mdm
        self.mX = inp.mX
        self.pol = inp.pol
        self.decay_mode = inp.decay_mode
        self.nsamples = inp.nsamples
        self.csv_path = inp.csv_path

    #Samples the appropriate 'box' shape of the final state spectrum. 
    def _sample_box(self, mdaughter2, cmcmm):
        if cmcmm == None:
            cmcmm = 0
        if self.XX:
            eta = arccosh(self.mdm / self.mX)
        else:
            eta = arccosh( (4 * self.mdm**2 + self.mX**2)/ (4 * self.mdm * self.mX))
        
        Emin = (self.mX**2 - mdaughter2**2) / (2 * self.mX) * (cosh(eta) - sinh(eta))
        Emax = (self.mX**2 - mdaughter2**2) / (2 * self.mX) * (cosh(eta) + sinh(eta))
        
        if self.pol == None:
            return self._sample_flat_box(Emin, Emax)
        elif self.pol in ['-0.5', '0.5']:
            return self._sample_ferm_box(Emin, Emax, cmcmm)
        elif self.pol in ['-1','0', '1', 'L', 'T']:
            return self._sample_vec_box(Emin, Emax, cmcmm)

    #Samples a flat box.
    def _sample_flat_box(self, Emin, Emax):
        u = rand()
        return Emin + u * (Emax - Emin)
        
    #Samples the box when X is polarized and fermionic.
    def _sample_ferm_box(self, Emin, Emax, cmcmm):
        if float(cmcmm) == 0:
            return self._sample_flat_box(Emin, Emax)
        if self.XX:
            Ex = self.mdm 
        else:
            Ex = self.mdm + self.mX**2 / (4 * self.mdm)
        return self._reject_sampling(Ex, 1 + cmcmm, self._ferm_func, cmcmm)

    #Returns the dNdE for the fermionic box.
    def _ferm_func(self, E, Ex, cmcmm):
        if self.pol == '-0.5':
            return 1 - cmcmm * (2*E - Ex)/(sqrt(Ex**2 - self.mX**2))
        else:
            return 1 + cmcmm * (2*E - Ex)/(sqrt(Ex**2 - self.mX**2))

    #Samples the box when X is polarized and vector-like.
    def _sample_vec_box(self, Emin, Emax, cmcmm):
        if self.XX:
            Ex = self.mdm 
        else:
            Ex = self.mdm + self.mX**2 / (4 * self.mdm)
        if self.pol in ['0', 'L']:
            ymax = self._vec_func(0.5 * Ex, Ex, cmcmm)
        else:
            ymax = max(self._vec_func(Emin, Ex, cmcmm), self._vec_func(Emax, Ex, cmcmm)) 
        return self._reject_sampling(Ex, ymax, self._vec_func, cmcmm)
    
    #Box shapes of a spin-1 X
    def _vec_func(self, E, Ex, cmcmm):
        if self.pol in ['0', 'L']:
            return -4*E**2 + 4*Ex*E - self.mX**2
        elif self.pol in ['-1','1']:
            return 4*E**2 - 4*Ex*E + 2*Ex**2 - self.mX**2 + sign(float(self.pol))*2*cmcmm*(2*E-Ex)*sqrt(Ex**2-self.mX**2)
        elif self.pol == 'T':
            return 4*E**2 - 4*Ex*E + 2*Ex**2 - self.mX**2

    #Performs rejections sampling
    def _reject_sampling(self, Ex, ymax, func, cmcmm):
        xmin = 0.5 * (Ex - sqrt(Ex**2 - self.mX**2))
        xmax = 0.5 * (Ex + sqrt(Ex**2 - self.mX**2))
        xx = uniform(xmin, xmax)
        yy = uniform(0, ymax)
        while yy > func(xx, Ex, cmcmm):
            xx = uniform(xmin, xmax)
            yy = uniform(0, ymax)
        return xx
    
    #Samples the peak of the final state spectrum.
    def _sample_peak(self):
        return (4 * self.mdm**2 - self.mX**2) / (4 * self.mdm)

    #Samples the final state spectrum of a given particle via its CDF
    def _sampler(self, CDF):
        u = uniform(0, CDF[1][-1])
        index = bisect(CDF[1], u)
        z = (u - CDF[1][index -1]) / (CDF[1][index] - CDF[1][index-1])
        E = CDF[0][index-1] + z * (CDF[0][index] - CDF[0][index-1])
        return E

    #Boosts a particle to the CM frame of the colliding DM particles.
    def _boost_particle(self, E_nu):
        if self.XX:
            eta = arccosh(self.mdm / self.mX)
        else:
            eta = arccosh((4 * self.mdm**2 + self.mX**2)/(4*self.mdm*self.mX))
        
        E1 = E_nu * (cosh(eta) - sinh(eta))
        E2 = E_nu * (cosh(eta) + sinh(eta))

        u = rand()

        return E1 + u * (E2 - E1)

    #Samples and boosts a final state particle from its spectrum
    def _sample_spect(self, CDF):
        E_r = self._sampler(CDF)
        E_b = self._boost_particle(E_r)
        return E_b

    #Counts the number of final state particles resulting from the decay of particle X
    def _number_of_final_state_particles(self):
        N = [0]
        for mode in self.decay_mode:
            N.append(N[-1] + mode.num1)
            N.append(N[-1] + mode.num2)
        N = N[1:]
        return N

    #Counts the number of final state particles per DM annihilation
    def final_state_particles_per_collision(self):
        N = self._number_of_final_state_particles()
        if self.XX:
            return 2*N[-1]
        else:
            return (N[-1] + 1)
    
    #Returns the normalization constant in order to obtain the spectrum per annihilation/decay
    def norm_const(self):
        return self.nsamples / self.final_state_particles_per_collision()

    #Samples the total spectrum. Note the 0.5 factor in front of the final_state_particles_per_collission when self.XX == True. This is because the sampling is done for a single X, but _number_of_final_state_particles gives the total number for 2 X.
    def _sample_total(self):
        number_of_final_state_particles = self._number_of_final_state_particles()
        if self.XX:
            u = uniform(0, 0.5*self.final_state_particles_per_collision())
            index = bisect(number_of_final_state_particles, u)
            #Index of which decay mode it is i.e. 1st, 2nd, 3rd, etc.
            i1 = int(floor(float(index) / 2.0))
            #Index of which particle of the decay mode it is, i.e. 1st or 2nd
            i2 = int(index%2)

            tmp_CDF = self.decay_mode[i1].CDF[i2]
            mdaughter1 = self.decay_mode[i1].get_mass(self.decay_mode[i1].daughter[i2])
            mdaughter2 = self.decay_mode[i1].get_mass(self.decay_mode[i1].daughter[((i2+1)%2)])
            cmcmm = self.decay_mode[i1].cmcmm

            if tmp_CDF == None:
                return self._sample_box(mdaughter2, cmcmm)
            else:
                return self._sample_spect(tmp_CDF)

        else:
            u = uniform(0, self.final_state_particles_per_collision())
            if u > number_of_final_state_particles[-1]:
                return self._sample_peak()
            else:
                index = bisect(number_of_final_state_particles, u)
                #Index of which decay mode it is i.e. 1st, 2nd, 3rd, etc.
                i1 = int(floor(float(index) / 2.0))
                #Index of which particle of the decay mode it is, i.e. 1st or 2nd
                i2 = int(index%2)

                tmp_CDF = self.decay_mode[i1].CDF[i2]
                mdaughter1 = self.decay_mode[i1].get_mass(self.decay_mode[i1].daughter[i2])
                mdaughter2 = self.decay_mode[i1].get_mass(self.decay_mode[i1].daughter[((i2+1)%2)])
                cmcmm = self.decay_mode[i1].cmcmm

                if tmp_CDF == None:
                    return self._sample_box(mdaughter2, cmcmm)
                else:
                    return self._sample_spect(tmp_CDF)

    #Samples the total final state spectrum Nsamples times.    
    def sample_spectrum(self):
        arr = []
        for _ in tqdm(range(self.nsamples)):
            arr.append(self._sample_total())
        if self.csv_path == 'plot':
            plt.hist(arr,200)
            plt.yscale('log')
            plt.xlabel(r'$E \quad [{\rm GeV}]$')
            plt.ylabel(r'$\frac{{\rm d} N_\nu}{{\rm d} E} \quad [{\rm GeV}^{-1}]$')
            if self.XX:
                if self.ann:
                    plt.title(r'${{\rm DM DM}} \to X X  \quad M_{{\rm DM}} = {0} \quad [{{\rm GeV}}] \quad M_X = {1} \quad [{{\rm GeV}}]$'.format(self.mdm, self.mX))
                else:
                    plt.title(r'${{\rm DM}} \to X X  \quad M_{{\rm DM}} = {0} \quad [{{\rm GeV}}] \quad M_X = {1} \quad [{{\rm GeV}}]$'.format(2*self.mdm, self.mX))
            else:
                if self.ann:
                    plt.title(r'${{\rm DM DM}} \to \nu X \quad M_{{\rm DM}} = {0} \quad [{{\rm GeV}}] \quad M_X = {1} \quad [{{\rm GeV}}]$'.format(self.mdm, self.mX))
                else:
                    plt.title(r'${{\rm DM DM}} \to \nu X \quad M_{{\rm DM}} = {0} \quad [{{\rm GeV}}] \quad M_X = {1} \quad [{{\rm GeV}}]$'.format(2*self.mdm, self.mX))
            plt.show()
        elif self.csv_path == 'logplot':
            plt.hist(log10(arr),200)
            plt.yscale('log')
            plt.xlabel(r'$\log_{10}(E) \quad [{\rm GeV}]$')
            plt.ylabel(r'$\frac{{\rm d} N_\nu}{{\rm d} E} \quad [{\rm GeV}^{-1}]$')
            if self.XX:
                if self.ann:
                    plt.title(r'${{\rm DM DM}} \to X X  \quad M_{{\rm DM}} = {0} \quad [{{\rm GeV}}] \quad M_X = {1} \quad [{{\rm GeV}}]$'.format(self.mdm, self.mX))
                else:
                    plt.title(r'${{\rm DM}} \to X X  \quad M_{{\rm DM}} = {0} \quad [{{\rm GeV}}] \quad M_X = {1} \quad [{{\rm GeV}}]$'.format(2*self.mdm, self.mX))
            else:
                if self.ann:
                    plt.title(r'${{\rm DM DM}} \to \nu X \quad M_{{\rm DM}} = {0} \quad [{{\rm GeV}}] \quad M_X = {1} \quad [{{\rm GeV}}]$'.format(self.mdm, self.mX))
                else:
                    plt.title(r'${{\rm DM DM}} \to \nu X \quad M_{{\rm DM}} = {0} \quad [{{\rm GeV}}] \quad M_X = {1} \quad [{{\rm GeV}}]$'.format(2*self.mdm, self.mX))
            plt.show()
        elif bool(self.csv_path):
            df = pd.DataFrame(arr)
            df.to_csv(f'{self.csv_path}', header = False, index = False)

#Contains all the relevant information for a decay mode of particle X
class Decay_modes():
    def __init__(self, frac, daughter1, daughter2, mX, final_state, cmcmm = None):
        self.frac = frac 
        self.daughter = [daughter1, daughter2]
        self.mX = mX
        self._check_decay_kinematics()
        self.final_state = final_state
        self.cmcmm = cmcmm
        self.CDF = [self._final_state_spectrum(daughter1, daughter2), self._final_state_spectrum(daughter2, daughter1)]
        self.num1 = self._final_state_counter(daughter1, daughter2)
        self.num2 = self._final_state_counter(daughter2, daughter1)

    #Gets the mass of the given daughter particle
    def get_mass(self, daughter):
        SM_dic = {'u': 0.002, 'd': 0.004, 's': 0.096, 'c': 1.28, 'b': 4.18, 't': 173.1, 'e': 0 ,'mu': 0.105, 'tau': 1.78, 'w': 80.4, 'z': 91.2, 'h':125.0, 'g': 0, 'ga': 0, 'v': 0}
        return SM_dic[daughter]
    
    #Checks if the decay of particle X is allowed
    def _check_decay_kinematics(self):
        if self.get_mass(self.daughter[0]) + self.get_mass(self.daughter[1]) > self.mX:
            raise Exception('The decay products of particle X are heavier than particle X itself; decay is thus kinematically forbidden. Aborting.')

    #Returns the cumulative distribution function of the final state spectrum of a given particle.
    def _final_state_spectrum(self, daughter1, daughter2):
        if daughter1 == 'v':
            return None
        elif daughter1 == 'ga' and self.final_state == 'Ga':
            md1 = self.get_mass(daughter1)
            md2 = self.get_mass(daughter2)
            mass = (self.mX**2 + md1**2 - md2**2) / (2 * self.mX)
            if mass < 5:
                raise Exception('One or more of the decay products of particle X has an equivalent mass lower than 5 GeV and can thus not be interpolated. Aborting.')
            data = Interpolate.annihilation_spectrum(mass, f'{daughter1}{daughter1}', self.final_state)
            CDF = [[],[]]

            num = 0
            x_old = 0.1 * min(data['x'])
            dndx_old = 0
            for x, dndx in zip(data['x'], data['dN/dx']):
                num += (dndx + 0.5 * abs(dndx - dndx_old)) * (1 - x_old / x)
                x_old = x
                dndx_old = dndx

                CDF[0].append(x * mass)
                CDF[1].append(0.5 * num / mass)

            CDF[0] = CDF[0][:-1]
            CDF[1] = CDF[1][:-1]
            
            CDF[0].append(1*mass)
            CDF[1].append(CDF[1][-1])

            CDF[0].append(mass)
            #The 1/mass factor is due to the rescaling of dN/dx to dN/dE
            CDF[1].append(CDF[1][-1]+2/mass)

            return CDF
        else:
            md1 = self.get_mass(daughter1)
            md2 = self.get_mass(daughter2)
            mass = (self.mX**2 + md1**2 - md2**2) / (2 * self.mX)
            if mass < 5:
                raise Exception('One or more of the decay products of particle X has an equivalent mass lower than 5 GeV and can thus not be interpolated. Aborting.')
            data = Interpolate.annihilation_spectrum(mass, f'{daughter1}{daughter1}', self.final_state)
            CDF = [[],[]]

            num = 0
            x_old = 0.1 * min(data['x'])
            dndx_old = 0
            for x, dndx in zip(data['x'], data['dN/dx']):
                num += (dndx + 0.5 * abs(dndx - dndx_old)) * (1 - x_old / x)
                x_old = x
                dndx_old = dndx

                CDF[0].append(x * mass)
                CDF[1].append(0.5 * num / mass)

            return CDF
    
    #Counts the number of final state particles coming from the fragmentation and hadronization of daughter 1
    def _final_state_counter(self, daughter1, daughter2):
        if daughter1 == 'v' and self.final_state in ['Nuel', 'Numu', 'Nuta']:
            return self.frac
        elif daughter1 == 'v' and self.final_state == 'Ga':
            return 0
        if daughter1 == 'ga' and self.final_state == 'Ga':
            return self.frac
        else:
            md1 = self.get_mass(daughter1)
            md2 = self.get_mass(daughter2)
            mass = (self.mX**2 + md1**2 - md2**2) / (2 * self.mX)
            data = Interpolate.annihilation_spectrum(mass, f'{daughter1}{daughter1}', self.final_state)
            num = 0
            x_old = 0.1 * min(data['x'])
            dndx_old = 0
            for x, dndx in zip(data['x'], data['dN/dx']):
                num += (dndx + 0.5 * abs(dndx - dndx_old))  * (1 - x_old / x)
                x_old = x
                dndx_old = dndx
            return num * self.frac * 1.5