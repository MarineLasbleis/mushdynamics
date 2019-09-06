# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.2.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import mushdynamics as md


# %% [markdown]
# # Mush dynamics and compaction
#
# ## How to use the package
#
#

# %% [markdown]
# ## Run one simulation

# %%
## The scripts here are also in the folder ./script/runs.py

def run(options):
    """ Run a model with the provided options"""
    print(options)
    model = md.growth.Compaction(md.mush.velocity_Sramek, **options)
    model.run()
    
def param_growth(r, exp, t_max, n=2, N_fig=5, basefolder="", R_init=1e-3, N_max=5000):
    dt = t_max/N_fig
    folder_name = basefolder+"/exp_{:.2e}_t_max_{:.2e}_radius_{:.2e}".format(exp, t_max, r)
    options = {'advection': "FLS",
                'n': n,
                'delta': 1.,
                'eta': 1.,
                'psi0': 1.,
                'psiN': 0.6,
                'phi_init': 0.4,
                'K0': 1.,
                'sign': 1.,
                'BC': "dVdz==0",
                'coordinates': "spherical",
                "growth_rate_exponent": exp,
                'filename': 'IC',
                'time_max': t_max,
                'dt_print': dt,
                'output': folder_name,
                "R_init": R_init*r,
                "N_init": max(5, int(N_max*R_init)),
                "Ric_adim": r}
    return options

options = param_growth(10., 1., 0.1, basefolder="./test/", R_init=5e-3, N_max=8000)
run(options)

# %% [markdown]
# output is then in the folder ./test/: 

# %%
from glob import glob
print(glob("./test/*"))
print(glob("./test/*/*"))

# %% [markdown]
# ## Some examples of figures from the package
#
#

# %% [markdown]
# ### Numerical scheme for the discretization of the equations

# %%
md.mush.schema()

# %% [markdown]
# ### Growth scenarios for the inner core

# %%
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # for the sake of lisibility, we remove the warnings here
    md.growth.plot_growth([6, 4])

# %%
