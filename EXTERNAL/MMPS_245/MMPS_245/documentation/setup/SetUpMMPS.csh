#!/bin/tcsh -f

##############################################################################################################
#   MMIL SetUpMMPS.csh file
#   Created   2010-02-15 dhagler
#   Last Mod: 2013-12-28 dhagler
##############################################################################################################

################################################################################
### Remove previous MMPS Setup
if ( `env | grep -c MMPSVER` == 1 ) then
#  echo "Removing previous MMPS Setup..."
  unsetenv MMPSVER
  unsetenv MMPS_PARENT
  unsetenv MMPS_DIR
  unsetenv MMPS_LIB
  unsetenv MMPS_PARMS
  unsetenv MMPS_MATLAB
  unsetenv MMPS_EXT
  unsetenv MMPS_EXTMAT
  unsetenv INTEL_LICENSE_FILE

  set path=`echo $path | tr ' ' '\n' | awk '$0 !~ "/MMPS"' | paste -sd' '`
  set path=`echo $path | tr ' ' '\n' | awk '$0 !~ "/packages/opt/intel/mkl"' | paste -sd' '`
  set path=`echo $path | tr ' ' '\n' | awk '$0 !~ "/packages/opt/intel/cce"' | paste -sd' '`
endif

################################################################################
### Supply the desired MMPS version string as an input. Use 'default' MMPS
### version of 236 if not supplied.
setenv MMPSVER 236

### Parse input
if ($#argv == 0 ) then
  echo "No arguments supplied... using default MMPS Version $MMPSVER"
else
  setenv MMPSVER $argv[1]
endif

################################################################################
### Default path additions

setenv MMPS_PARENT ${PUBSW}/packages/MMPS

setenv MMPS_DIR ${MMPS_PARENT}/MMPS_${MMPSVER}
setenv MMPS_EXT ${MMPS_DIR}/external
setenv MMPS_LIB ${MMPS_DIR}/lib
setenv MMPS_PARMS ${MMPS_DIR}/parms
setenv MMPS_MATLAB ${MMPS_DIR}/matlab
setenv MMPS_EXTMAT ${MMPS_EXT}/matlab

set addpathlist = ( \
  ${MMPS_DIR}/bin \
  ${MMPS_DIR}/csh \
  ${MMPS_DIR}/ksh \
  ${MMPS_DIR}/perl \
  ${MMPS_DIR}/sh \
  ${MMPS_DIR}/tcl \
  $PUBSW/packages/dcmtk/3.6.0/bin \
  ${MMPS_EXT}/OpenMEEG/bin \
)
foreach dir ( $addpathlist )
  if (-e $dir && "$path" !~ *$dir*) set path = ($dir $path)
end
unset addpathlist dir

#########################################################################

# libraries required for Dominic Holland's unwarpB0, reg, etc.
if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH   ${MMPS_LIB}/BasicStructs:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH   ${MMPS_LIB}/BasicStructs
endif
setenv LD_LIBRARY_PATH   ${MMPS_LIB}/Interpolation:${LD_LIBRARY_PATH}
set mkl_csh = ${PUBSW}'/packages/opt/intel/mkl/10.0.2.018/tools/environment/mklvarsem64t.csh'
set iccenv_csh = ${PUBSW}'/packages/opt/intel/cce/10.1.013/bin/iccvars.csh'
if (-e $mkl_csh && !($?INTEL_LICENSE_FILE)) then
  source $mkl_csh
endif
if (-e $iccenv_csh && !($?INTEL_LICENSE_FILE)) then
  source $iccenv_csh
endif
if (-e ${MMPS_LIB}'/gsl') then
  setenv LD_LIBRARY_PATH   ${MMPS_LIB}/gsl/lib:${LD_LIBRARY_PATH}
endif

# for OpenMEEG
set omlib = ${MMPS_EXT}/OpenMEEG/lib
if (-e $omlib) then
  setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH":"$omlib
endif

################################################################################
### Clean path
set path=`echo $path | awk '{for(i=1;i<=NF;i++){if(!($i in a)){a[$i];printf s$i;s=" "}}}'`
