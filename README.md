# Dynamics of mush


## Overview


## Requirements and installation

We used Anaconda to package and install the package, but you can also directly install python3 and the packages through pip3. 

Required packages are: numpy scipy pandas matplotlib pyyaml

To be able to run the notebooks, please install jupyter notebook

For a complete installation through anaconda, you can create a new virtual environment with the required packages:

```
conda create -n name_env python=3.6

conda activate name_env

conda install numpy scipy pandas matplotlib pyyaml jupyter
```
And finally, to install mushdynamics, while in the folder:

```
python setup.py install
```

To be able to use git with the jupyter notebook, we use [jupytext](https://github.com/mwouts/jupytext).

## Tests

There is no unit test implemented yet. But there are one script with a list of small programs to visually test if the solvers are doing correct things and to compare the advection solvers. See tests.py . 