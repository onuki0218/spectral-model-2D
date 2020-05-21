#!/bin/sh

#------------- pjsub option ---------------
#PJM -L rscgrp=debug-cache
#PJM -L node=1
#PJM --mpi proc=5
#PJM --omp thread=4
#PJM -L elapse=0:30:00
#PJM -g gi55
#PJM -j
#PJM -s
#PJM -m b,e,r
#PJM --restart

#----------- Pinning setting --------------
source /usr/local/bin/mpi_core_setting.sh

#---------- Program execution -------------
data_dir='../data/'
output_dir=${data_dir}data${PJM_JOBID}/
mkdir ${output_dir}
echo Output directory: ${output_dir}
export setting="setting"
export dir_name_result=${output_dir}
module load fftw
module load mpi-fftw
mpiexec.hydra -n ${PJM_MPI_PROC} ./main.exe
