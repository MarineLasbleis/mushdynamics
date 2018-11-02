""" From Huguet et al. 2018

Calculation of the growth rate and radius of the IC with supercooling.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.optimize import fsolve

year = 365.25*3600*24


class Evolution:
    """ Evolution of inner core.

    Return the radius and growth rate as function of time.
    """

    def __init__(self, delta_T, Qcmb, N_time=1e5):
        self.Qcmb = Qcmb
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
        self.N_time = N_time
        self.file = "supercooling_{:.2e}_{:.2e}.txt".format(delta_T, Qcmb)

    def initilization(self):
        self.time = np.zeros(int(self.N_time))
        self.dt = np.zeros_like(self.time)
        self.dric = np.zeros_like(self.time)
        self.ric = np.zeros_like(self.time)
        self.Tic = np.zeros_like(self.time)
        self.supercooling = np.zeros_like(self.time)
        self.Tis_center = np.zeros_like(self.time)
        self.Tcmb = np.zeros_like(self.time)

    def new_dt(self, dot_r, dot_T, dr0=1.):
        # dr0 = 50 # max variation of r allowed
        # dT0 = 0.1 # max variation of Tcmb allowed
        dt0_r = np.abs(dr0 /dot_r) *0.5
        return min(1000*year, dt0_r)

    def run_constant_temperature(self):
        self.initilization()
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

    def run(self):
        self.initilization()
        T_center = self.T_Fe(0., 0.) - self.delta_T # initial temperature at center.
        self.dric[0] = 0. #self.G*(self.delta_T)**2
        self.Tic[0] = T_center
        self.supercooling = self.delta_T
        self.Tis_center[0] = T_center
        self.Tcmb[0] = self.T_is(self.rc, T_center)
        self.dt[0] = 1e-2*year
        for i, t in enumerate(self.time[1:]):
            dric, Tic = self.icoregrowth(self.Tis_center[i], self.ric[i], 0.)
            self.dric[i+1] = dric
            self.Tic[i+1] = Tic
            self.time[i+1] = self.time[i] + self.dt[i]
            self.dt[i+1] = self.new_dt(dric, 0.)
            self.ric[i+1] = self.ric[i] + dric*self.dt[i]
            #self.supercooling[i] = np.abs(self.T_Fe(self.ric[i]) - Tic)
            self.Tis_center[i+1] = self.Tis_center[i] + self.dTcmb(dric, self.ric[i])*self.dt[i]
            self.Tcmb[i+1] = self.T_is(self.rc, self.Tis_center[i+1]) # TODO change to r_outer_core
            if i%1e4==0: print("Iteration: {}, Time: {:.2e}year, dt: {:.2e}year,  size IC: {:.2f}km, dr: {:.3f}m.".format(i, self.time[i+1]/year, self.dt[i]/year,  self.ric[i+1]/1e3, dric))

    def run_print(self):
        """ Run and print directly in a text file """
        # initialisation
        T_center = self.T_Fe(0., 0.) - self.delta_T
        Tic = T_center
        time = 0.
        Tis_center = T_center
        ric = 0.
        dric = 0.
        i = 0
        dt = 1e-2*year
        n_i = 1e4

        with open(self.file, 'w') as f:
            f.write("iteration time(year) dt(year) ric(m) dric(m) \dot{r}_ic(m/s) T_center(K)\n")
            f.write("{} {:6.3e} {:6.3e} {:6.3f} {:6.3e} {:6.3e} {:6.2f}\n".format(i, time/year, dt/year, ric, dric, dric/dt,Tis_center))

        while i<self.N_time and ric<1221e3:
            dric, Tic = self.icoregrowth(Tis_center, ric, 0.)
            Tis_center += self.dTcmb(dric, ric)*dt
            ric += dric*dt
            time += dt
            i +=1
            dt = self.new_dt(dric, 0.)
            if i%1e5==0: print("Iteration: {}, Time: {:.2e}year, dt: {:.2e}year,  size IC: {:.2f}km, dr: {:.3e}m.".format(i, time/year, dt/year,  ric/1e3, dric))
            if i%n_i==0:
                with open(self.file, 'a') as f:
                    f.write("{} {:6.3e} {:6.3e} {:6.3f} {:6.3e} {:6.3e} {:6.2f}\n".format(i, time/year, dt/year, ric, dric, dric/dt, Tis_center))


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
        a, _, ierr, msg = fsolve(func, np.array([1e-9, self.T_is(ric, T)]), xtol=1e-9, args=(T, ric, cc,), full_output=True)
        if ierr != 1:
            print(ric, ierr, a, msg)
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
        ax[0].plot(self.Tic[:], self.ric[:], 'k', label="T_ic")
        ax[0].plot(self.Tis_center[:], self.ric[:], 'g', label="Tcenter")
        ax[2].plot(self.dric[:], self.ric[:], 'g', label="dot r")
        ax[0].legend()
        ax[1].plot(self.time[:]/year, self.ric[:])
        ax[0].set_ylabel("Radius")
        ax[0].set_xlabel("Temperature")
        ax[1].set_xlabel("Time")
        ax[1].set_xscale("log")
        print("Time: {:.2e}year, size IC: {:.2f}m.".format(self.time[-1]/year, self.ric[-1]))


def load_data(file):
    data = pd.read_csv(file, sep=" ", index_col=False)
    return data

def fig_from_file(file):
    data = load_data(file)
    fig, ax = plt.subplots()
    data.plot(x="time(year)", y="ric(m)", ax=ax, logx=True)
    return data

if __name__ == "__main__":

    delta_T = 30 #initial supercooling

    test = Evolution(100, 1e13, 5e6)
    test.run_print()
    fig_from_file(test.file)
    test = Evolution(30, 1e13, 3e6)
    test.run_print()
    fig_from_file(test.file)
    test = Evolution(120, 1e13, 3e6)
    test.run_print()
    fig_from_file(test.file)
    test = Evolution(1, 1e13, 5e6)
    test.run_print()
    fig_from_file(test.file)
    test = Evolution(100, 1e14, 5e6)
    test.run_print()
    fig_from_file(test.file)

    #test.run_print()
    # test.plot()

    #test = Evolution(delta_T, 1e9)
    #test.run()
    #test.plot()
    #test.run_constant_temperature()
    #test.plot()

    

    plt.show()
