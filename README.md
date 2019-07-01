# TAGSS

Supplementary package for the summer school [Trieste Algebraic Geometry Summer School (TAGSS) 2019 - Algebraic Geometry towards Applications](http://indico.ictp.it/event/8695/).

To get things started clone or download this repository and switch to this folder and instantiate the project.

```
# clone repository
git clone https://github.com/saschatimme/TAGSS.git
# switch to the folder
cd TAGSS
# install all packages necessary
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
# start the repository and open the provided Jupyter notebook
julia --project -e 'using TAGSS; ed_notebook();'
```
