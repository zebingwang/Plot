#!/bin/bash
/bin/hostname
gcc -v
pwd
export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/
source /cvmfs/cms.cern.ch/cmsset_default.sh


#### job

##########hep_sub runjob.sh -g cms -mem 8000 -wt mid -o job.out -e job.err
#python ALP_plot_param.py -m -y 2016 

#python ALP_Optimization.py -y run2 -o ./optimize_run2 --doOpt -c 5
#python ALP_plot_param.py -y run2 -m --ln

#python ALP_plot_param.py -y run2 -m -S #--ln #--cut --mA M30
#python ALP_plot_param.py -y run2 -m 
python ALP_plot_bkgCorr.py -y run2 -m