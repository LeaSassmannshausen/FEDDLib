#!/bin/bash -l
#
#SBATCH --job-name=sci_test
#SBATCH --comment="sci_test"
#SBATCH --time=24:00:00         ### time the job will appr. run
#SBATCH --nodes=4               ### Node count required for the job
#SBATCH --output=/home/hpc/k105be/k105be13/slurm/outputs/sci_test.%j.out              ### output file for console output
#SBATCH --error=/home/hpc/k105be/k105be13/slurm/outputs/sci_test.%j.err               ### output file for console error
#SBATCH --ntasks=150              ### Number of tasks per job    (or next line, should one of them)

unset SLURM_EXPORT_ENV

#source $HOME/opt/spack/share/spack/setup-env.sh
module load intel/2023.2.1 intelmpi/2021.10.0 mkl/2023.2.0 
#hdf5/1.12.2-oneapi2023.2.0-impi-6zjqfkp metis/5.1.0-oneapi2023.2.0-mrb7o76 parmetis/4.0.3-oneapi2023.2.0-impi-yofi5yp boost/1.79.0-oneapi2023.2.0-impi-sucm53o ninja/1.11.1-gcc8.5.0-ckyge76 cmake
#export I_MPI_PMI_LIBRARY=/lib64/libpmi2.so

#./do-config-trilinos.sh
#ninja
#ctest -j 30
#ninja install
#srun --mpi=mpi2 ./problems_sci_Test.exe
srun ./problems_sci_Test.exe