""" Series of runs """

print("bouh")
import mush
print("mush imported")
import growth
print("growth imported")

if __name__ == "__main__":

    r_max = [200.]# [0.2, 0.01, 0.1]# , 2., 5., 10., 20., 50., 100.]
    N_fig = 10
    exp_velocity = [0.5] # , 0.5, 1./3.]
    coeff_velocity = [10., 1., 0.1 ]#[1e-1, 1e-2, 1e-3, 1e-4, 1., 10., 20., 100.] # [1., 0.5, 0.05]# [0.05, 0.1, 1., 2., 5., 10., 20., 50., 100.]
    r0_supercooling = [0.65, 0.6]

    def new_options(**param):
        r_max = 10.
        t_max = (10/2.)**2
        dt = t_max/20.
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
                'filename': 'IC_Sramek_cart',
                'time_max': t_max,
                'dt_print': dt,
                'coeff_velocity': 2.,
                "R_init": 0.01,
                "N_init": max(5, int(2000./r_max))}
        options = {**options, **param}
        return options

    for r in r_max:
        for exp in exp_velocity:
            for coeff in coeff_velocity:
                t_max = (r/coeff)**(1/exp)
                print(t_max)
                dt = t_max/N_fig
                print(dt)        
# r0_supercooling = 0.8*r
# folder_name = "compaction/exp_{:.2f}_coeff_{:.2f}_radius_{:.2f}".format(exp, coeff, r)
                folder_name = "compaction/exp_{:.2f}_coeff_{:.4f}_radius_{:.2f}".format(exp, coeff, r )
                options = new_options(growth_rate_exponent=exp,
                                        time_max=t_max,
                                        dt_print=dt,
                                        output=folder_name,
                                        coeff_velocity=coeff, 
                                        R_init = 0.005*r)
                print(folder_name)
                print("Time to be computed: {:.2e}, dt for print: {:.2e}".format(t_max, dt))
                #Model = growth.Compaction(mush.velocity_Sramek, **options)
                #options["output"] = "compaction_supercooling/exp_{:.2f}_coeff_{:.2f}_radius_{:.2f}".format(exp, coeff, r)
                #options["t0_supercooling"] = 1.
                #options["Dt_supercooling"] = 60.
                #options["output"] = "compaction_supercooling/exp_{:.2f}_coeff_{:.2f}_radius_{:.2f}_Dt_{:.2f}".format(exp, coeff, r, options["Dt_supercooling"])
                #Model = growth.Compaction_Supercooling(mush.velocity_Sramek, **options)
                Model = growth.Compaction(mush.velocity_Sramek, **options)
                Model.run()
