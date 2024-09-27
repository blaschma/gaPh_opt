# Phonon Transport
Calculates the phononic transport according to [[1]](#1).
## Requirements
* Python 3.x
* [tmoutproc](https://github.com/blaschma/tmoutproc/)


## Calculation of phonon transmission
```` 
python3 phonon_transport.py config_file
```` 
### Preparation
* Build molecule in Avogadro. Save as .xyz file 
* Relaxation and hessian:
    * Geometry optimization e.g. using xtb (https://xtb-docs.readthedocs.io/en/latest/optimization.html)
    * Align molecule -> x axis through anchoring atoms (not hydrogen). Important if in_plane option is used
    * Calculate hessian e.g using xtb (https://xtb-docs.readthedocs.io/en/latest/hessian.html)

### Config file
```` 
[Data Input]
data_path= #path where data is located
hessian_name=hessian #filename of hessian
coord_name=coord.xyz #filename of coord file (.xyz and Turbomole format allowed)

[Calculation]
n_l=5,1,2 #indices of atoms (-> ordering of coord file) connected to left lead
n_r=1,2,3 #indices of atoms (-> ordering of coord file) connected to right lead
M_L=Au # atom type in left/right lead
M_C=Au # atom type in center connected to lead
gamma= #coupling constant from [Phonon interference effects in molecular junctions](https://doi.org/10.1063/1.4849178)
E_D= #Debey energy
N= #number of grid points for transmission calculation
in_plane=False #in_plane option from [Phonon interference effects in molecular junctions](https://doi.org/10.1063/1.4849178)
T_min= #lower bound for thermal conductance integral (avoid zero)
T_max= #upper bound for thermal conductance integral
kappa_grid_points= #number of grid point in thermal conductance integral

[Eigenchannel]
eigenchannel=True (True: Eigenchannels are calculated. See comment)
every_nth=1 (for all, -1 for none. Specifies number of plotted eigenchannels in data_path/transport_channels.pdf)
channel_max=3 (number of plotted and stored eigenchannels. All channels are calculated)

[Data Output]
plot_g=True #plot surface green function 

````

### Output
* data_path/phonon_trans.dat
* data_path/kappa.dat
* data_path/transport.pdf
* data_path/g0.pdf (optional, see plot_g)
* data_path/transport_channels.pdf (optional, see Eigenchannel)
* data_path/transport_channels.dat (NOT IMPLEMENTED YET, see Eigenchannel)
* data_path/eigenchannels/*.nmd (optional, see Eigenchannel)

> **_NOTE:_**  transport_channels.dat not implemented yet.

### Calculation of Eigenchannels
If eigenchannel=True is set, the transmission eigenchannels according to  [[2]](#1) are calculated. The total transmission is then calculated as sum over all eigenchannels. The first channel_max channels are plotted in data_path/transport_channels.pdf. every_nth specifies which eigenchannels are written to *nmd file. 

### Calculation of Propagator elements
> **_NOTE:_**  Not implemented yet.

### Calculation of Participation Ratio
> **_NOTE:_**  Not implemented yet.

## Calculation of thermal conductance (Standalone)
### Usage
```` 
python3 calculate_kappa.py config_file
```` 
Calculates thermal conductance from phonon transmission. Energy must be in Hartrees!
### Preperation
Transport calculation for transmission (does not necessarily have to be calculated with this program)
### Config file
A reduced config file is sufficient for this
```` 
[Data Input]
data_path= #path where data is located
transp_name= #name of file containing phonon transmission
transp_units = [har],[sqrt(har/(bohr**2*u))] # hartree is default

[Calculation]
kappa_int_lower_E=0 #lower integral limit in kappa Energy integral in meV (optional, for further analysis). See commet below
kappa_int_upper_E=5 #upper integral limit in kappa Energy integral in meV (optional, for further analysis). See comment below
T_min= #lower bound for thermal conductance integral (avoid zero)
T_max= #upper bound for thermal conductance integral
kappa_grid_points= #number of grid point in thermal conductance integral
````
If kappa_int_lower_E and kappa_int_upper_E are set, the cumulative thermal conductance $\kappa^{\mathrm{c}}_{\mathrm{ph}}$ is calculated with this integral limits at 300K.

### Output
* data_path/kappa.dat
* data_path/kappa.pdf
* data_path/kappa_c.dat (optional)
* data_path/kappa_c.pdf (optional)

## Calculation of Phonon Eigenchannels (Standalone)
> **_NOTE:_**  Not implemented yet.


## Planned features
* Eigenchannel
  * .g98 files for Eigenchannels
  * Calculation of Phonon Eigenchannels (Standalone)
  * Writeout of eigenchannel data
  
* Multiple electrode models
* Database for coupling parameters or more consistent calculation
* Example files
* Writeout of propagator elements
* Writeout of Participation ratio



## References
<a id="1">[1]</a> 
Markussen, T. (2013).
Phonon interference effects in molecular junctions. 
The Journal of chemical physics, 139(24), 244101.
[https://doi.org/10.1063/1.4849178](https://doi.org/10.1063/1.4849178) \
<a id="1">[2]</a> 
Klöckner, J. C., Cuevas, J. C., & Pauly, F. (2018). Transmission eigenchannels for coherent phonon transport. Physical Review B, 97(15), 155432.
[https://doi.org/10.1103/PhysRevB.97.155432](https://doi.org/10.1103/PhysRevB.97.155432)

***
Matthias Blaschke [matthias.blaschke@physik.uni-augsburg.de](matthias.blaschke@pyhsik.uni-augsburg.de)
