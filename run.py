""" Series of runs """


import mush
import growth

if __name__ == "__main__":

    r_max = [200.]# [0.2, 0.01, 0.1]# , 2., 5., 10., 20., 50., 100.]
    N_fig = 10
    exp_velocity = [0.5] # , 0.5, 1./3.]
    coeff_velocity = [10., 1., 0.1 ]#[1e-1, 1e-2, 1e-3, 1e-4, 1., 10., 20., 100.] # [1., 0.5, 0.05]# [0.05, 0.1, 1., 2., 5., 10., 20., 50., 100.]
    r0_supercooling = [0.65, 0.6]

    def new_options(**param):
        r_max = 10.
        t_max = (10/2.)**2
        dt = t_max/20
        options = {'advection': "FLS",
                'delta': 1.,
                'eta': 1.,
                'psi0': 1.,
                'psiN': 0.6,
                'phi_init': 0.4,
                'K0': 1.,
                'delta': 1.,
                'sign': 1.,
                'BC': "dVdz==0",
                'coordinates': "spherical",
                "growth_rate_exponent": 0.5,
                'filename': 'IC_Sramek',
                'time_max': t_max,
                'dt_print': dt,
                'coeff_velocity': 2.,
                'output': "compaction/",
                "R_init": 0.001,
                "N_init": 5} # min(5, int(4e3/r_max))}
        options = {**options, **param}
        return options

    for r in r_max:
        for exp in exp_velocity:
            for r0 in r0_supercooling: 
                for coeff in coeff_velocity:
                    t_max = (r/coeff)**(1/exp)
                    dt = t_max/N_fig
                    folder_name = "data_for_fig_05/exp_{:.2f}_coeff_{:.5f}_radius_{:.2f}".format(exp, coeff, r)
                    options = new_options(growth_rate_exponent=exp,
                                            time_max=t_max,
                                            dt_print=dt,
                                            output=folder_name,
                                            coeff_velocity=coeff, 
                                            R_init = 0.0005*r)
                    print(folder_name)
                    print("Time to be computed: {:.2e}, dt for print: {:.2e}".format(t_max, dt))
                    #Model = growth.Compaction(mush.velocity_Sramek, **options)
                    # options["output"] = "compaction_supercooling/exp_{:.2f}_coeff_{:.2f}_radius_{:.2f}".format(exp, coeff, r)
                    options["t0_supercooling"] = 0.01*t_max
                    options["r0_supercooling"] = r0*r
                    options["output"] = "data_for_fig_supercooling/exp_{:.2f}_coeff_{:.2f}_radius_{:.2f}_r0_{:.2f}".format(exp, coeff, r, options["r0_supercooling"])
                    Model = growth.Compaction_Supercooling(mush.velocity_Sramek, **options)
                    #Model = growth.Compaction(mush.velocity_Sramek, **options)
                    Model.run()
