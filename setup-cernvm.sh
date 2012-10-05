
# Setup root and python with gcc 43 because the build fails
# with gcc46

#source /cvmfs/sft.cern.ch/lcg/external/gcc/4.3.5/x86_64-slc5-gcc34-opt/setup.sh
#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.00/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh 
#PYTHONHOME=/cvmfs/sft.cern.ch/lcg/external/Python/2.6.5/x86_64-slc5-gcc43-opt
#export PATH=${PYTHONHOME}/bin:${PATH}
#export LD_LIBRARY_PATH=${PYTHONHOME}/lib:${LD_LIBRARY_PATH}

# It turns out gcc 46 is ok after all
source /cvmfs/sft.cern.ch/lcg/hepsoft/0.6/x86_64-slc5-gcc46-opt/setup.sh

