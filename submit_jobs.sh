#!/bin/bash

problem_sizes=(1200 2400 3600 4800 6000 7200 8400 9600)

rank_counts=(2 4 5 6 8 10 12 15 20 30 40 60)

for size in "${problem_sizes[@]}"; do
    for ranks in "${rank_counts[@]}"; do
        # Each node has 20 physical cores (2 CPUs * 10 cores each)
        nodes_needed=$(( (ranks + 19) / 20 )) # Ensure at least 1 node is used

        echo "Problem size: $size, Rank count: $ranks, Nodes needed: $nodes_needed"

        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=heat_eq_${size}_${ranks}
#SBATCH --output=heat_eq_%j.out
#SBATCH --partition=peda
#SBATCH --nodes=$nodes_needed
#SBATCH --ntasks-per-node=20 # Max out the cores per node first
#SBATCH --time=01:00:00

#module load mpi
mpirun -np $ranks /scratch/pipebay/f4-hpc/isotherm/build/heat2d_parallel $size
EOF
    done
done

