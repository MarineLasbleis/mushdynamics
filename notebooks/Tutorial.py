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
# ## Some examples of figures from the package
#
# ### Growth scenarios for the inner core

# %%
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # for the sake of lisibility, we remove the warnings here
    md.growth.plot_growth([10, 8])

# %%
