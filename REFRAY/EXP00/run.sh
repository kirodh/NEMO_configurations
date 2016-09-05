## ======================================================================
## Path ... should be completed
## ======================================================================
MAINPATH=$HOME             # path of NEMO trunk

DIR_ARCHI_INIT=$MAINPATH/NEMOGCM/CONFIG/REFRAY/EXP00/input     # PATH for the initial condition file
DIR_ARCHI_FORC=$MAINPATH/NEMOGCM/CONFIG/REFRAY/EXP00/input     # PATH for the forcing input files


## ======================================================================
## =======  Nothing to be changed bellow this line  =====================
## ======================================================================

mkdir rundir
cd rundir
## ======================================================================
## file needed for the run
## ======================================================================
cp ${MAINPATH}/NEMOGCM/CONFIG/SHARED/namelist_ref . 
cp ../*.xml            .
cp ../namelist_cfg     .
cp ../opa              .

## ======================================================================
## Climatology of temperature and chlorophyll
## ======================================================================
ln -s ${DIR_ARCHI_INIT}/init_PAPASTATION_m06d15.nc  ./init_PAPASTATION_m06d15.nc
ln -s ${DIR_ARCHI_INIT}/chlorophyll_PAPASTATION.nc  ./chlorophyll_PAPASTATION.nc

## ======================================================================
## Atmospherical forcings
## ======================================================================
ln -s ${DIR_ARCHI_FORC}/forcin* .

./opa






