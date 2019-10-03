import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mushdynamics",
    version="0.1",
    author="Marine Lasbleis",
    author_email="marine.lasbleis@gmail.com",
    description="Solvers for 2-phase flow dynamics of a compacting mush",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marinelasbleis/mushdynamics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

