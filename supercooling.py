""" From Huguet et al. 2018

Calculation of the growth rate and radius of the IC with supercooling.
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve

year = 365.25*3600*24


class Evolution:
    """ Evolution of inner core.

    Return the radius and growth rate as function of time.
    """

    def __init__(self, delta_T, Qcmb):
        self.Qcmb = Qcmb
        print(self.Qcmb)
        self.delta_T = delta_T
        # parameters
        self.longueur = 6500e3 # Lp = sqrt(3Cp/2pialpharhoG)
        self.rc = 3600e3 # radius of CMB
        self.G = 1e-4
        self.rhol = 12500
        self.rho0 = 11900
        self.u = 3e-4
        self.cp = 785
        self.LFe = 750e3 # latent heat of crystallization
        self.TFe0 = 6500 #melting T of iron at center
        self.gamma = 1.5 # Gruneisen parameter
        self.eta_c = 0.79 # parameters (see after eq. A.3 in Huguet et al 2018)
        self.Mc = 1.95e24
        self.Eg =  3e5
        # numerical solver
        #time_max = 1e6
        N_time = 1e6
        self.time = np.zeros(int(N_time)) # np.linspace(0, time_max*year, N_time)
        # self.time_years = self.time/year
        self.dt = np.zeros_like(self.time)
        # self.dt = self.time[1] - self.time[0] # 1ks #to be determined automatically
        self.dric = np.zeros_like(self.time)
        self.ric = np.zeros_like(self.time)
        self.Tic = np.zeros_like(self.time)
        self.supercooling = np.zeros_like(self.time)
        self.Tis_center = np.zeros_like(self.time)
        self.Tcmb = np.zeros_like(self.time)

    def reinitilization(self):
        self.dric = np.zeros_like(self.time)
        self.ric = np.zeros_like(self.time)
        self.Tic = np.zeros_like(self.time)
        self.supercooling = np.zeros_like(self.time)
        self.Tis_center = np.zeros_like(self.time)
        self.Tcmb = np.zeros_like(self.time)

    def new_dt(self, dot_r, dot_T):
        dr0 = 10 # max variation of r allowed
        # dT0 = 0.1 # max variation of Tcmb allowed 
        dt0_r = np.abs(dr0 /dot_r) *0.5
        return min(self.dt[0], dt0_r)

    def run_constant_temperature(self):
        self.reinitilization()
        T_center = self.T_Fe(0., 0.) - self.delta_T # initial temperature at center. 
        self.dric[0] = 0. #self.G*(self.delta_T)**2
        self.Tic[0] = T_center
        self.supercooling = self.delta_T
        self.Tis_center[:] = T_center
        self.Tcmb[0] = self.T_is(self.rc, T_center)
        self.dt[0] = 0.1*year

        for i, t in enumerate(self.time[:-1]):
            dric, Tic = self.icoregrowth(T_center, self.ric[i], 0.)
            self.dric[i+1] = dric
            self.Tic[i+1] = Tic
            self.ric[i+1] = self.ric[i] + dric*self.dt[i]
            self.time[i+1] = self.time[i] + self.dt[i]
            self.dt[i+1] = self.new_dt(dric, 0.)
            if i%1e4==0: print(i, self.time[i+1]/year)

            #self.supercooling[i] = np.abs(self.T_Fe(self.ric[i]) - Tic)

    def run(self):
        self.reinitilization()
        T_center = self.T_Fe(0., 0.) - self.delta_T # initial temperature at center.
        self.dric[0] = 0. #self.G*(self.delta_T)**2
        self.Tic[0] = T_center
        self.supercooling = self.delta_T
        self.Tis_center[0] = T_center
        self.Tcmb[0] = self.T_is(self.rc, T_center)
        self.dt[0] = 10*year
        for i, t in enumerate(self.time[1:]):
            dric, Tic = self.icoregrowth(self.Tis_center[i], self.ric[i], 0.)
            self.dric[i+1] = dric
            self.Tic[i+1] = Tic
            self.time[i+1] = self.time[i] + self.dt[i]
            self.dt[i+1] = self.new_dt(dric, 0.)
            self.ric[i+1] = self.ric[i] + dric*self.dt[i]
            #self.supercooling[i] = np.abs(self.T_Fe(self.ric[i]) - Tic)
            self.Tis_center[i+1] = self.Tis_center[i] + self.dTcmb(dric, self.ric[i+1])*self.dt[i]
            self.Tcmb[i+1] = self.T_is(self.ric[i+1], self.Tis_center[i+1])
            if i%1e4==0: print(i, self.time[i+1]/year)

    def T_is(self, radius, T0):
        # T0 is temperature at center
        return T0*np.exp(-(radius**2)/self.longueur**2)

    def T_Fe(self,radius, xi=0.):
        return self.TFe0*np.exp(-2*(1-1/3/self.gamma)*radius**2/self.longueur**2) - xi

    def icoregrowth(self, T, ric, cc):
        # Eq. A.5 and A.6 in Huguet+2018
        # returns dric/dt in m/s
        def func(y, T, ric, cc):
            rd = y[0]
            Ti = y[1]
            Ta = self.T_is(ric, T)
            Tm = self.T_Fe(ric, cc)
            y0 = rd - self.G*(Tm - Ti)**2
            y1 = self.rho0*(self.cp*(Tm - Ta) + self.LFe)*rd - self.rhol*self.cp*self.u*(Ti - Ta)
            return np.array([y0, y1])
        a, _, ierr, msg = fsolve(func, np.array([1e-9, T]), xtol=1e-4, args=(T, ric, cc,), full_output=True)
        if ierr != 1:
            print(ierr, a, msg)
        return a[0], a[1]

    def QLQG(self, dricdt, ric):
        """ term Q_L + Q_G with equation A.5 Huguet (2018) """
        Aic = 4.*np.pi*ric**2
        return (self.LFe + self.Eg)*self.rho0*dricdt*Aic

    def dTcmb(self, dricdt, ric):
        return self.eta_c/self.Mc/self.cp * (self.QLQG(dricdt, ric) - self.Qcmb)

    def plot(self):
        fig, ax = plt.subplots(1,3, sharey=True)
        r = np.linspace(0, 1221e3, 30)
        ax[0].plot(self.T_is(r, 6400), r, 'r', label="T isentrope_beg")
        #ax[0].plot(self.T_is(r, self.Tis_center[-1]), r, 'r', label="T isentrope_end")
        #ax[0].plot(self.T_is(r, self.Tis_center[500]), r, 'r', label="T isentrope_middle")
        ax[0].plot(self.T_Fe(r), r, 'b', label="melting T of iron")
        ax[0].plot(self.Tic[:-1], self.ric[:-1], 'k', label="T_ic")
        ax[0].plot(self.Tis_center[:-1], self.ric[:-1], 'g', label="Tcenter")
        ax[2].plot(np.diff(self.Tcmb)[1:], self.ric[1:-1], 'g', label="dot r")
        ax[0].legend()
        ax[1].plot(self.time[:-1]/year, self.ric[:-1])
        ax[0].set_ylabel("Radius")
        ax[0].set_xlabel("Temperature")
        ax[1].set_xlabel("Time")
        ax[1].set_xscale("log")
        print(self.time[-1]/year, self.ric[-1])


if __name__ == "__main__":

    delta_T = 100 #initial supercooling

    test = Evolution(delta_T, 1e13)
    test.run()
    test.plot()

    #test = Evolution(delta_T, 1e9)
    #test.run()
    #test.plot()
    #test.run_constant_temperature()
    #test.plot()
    plt.show()
