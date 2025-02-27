#!/bin/bash
#SBATCH --nodes 1
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J PARAM_SPACE
#SBATCH -t 0:40:0
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=siennaw@berkeley.edu

##### SBATCH -account m1266

# Define the input CSV file
CSV_FILE="parameter_space_ws_pmax.csv"

# Skip the header and read the CSV line by line
 tail -n +1 "$CSV_FILE" | while IFS="," read -r ws1 ws2 pmax1 pmax2 fout_name; do
    # Call the Python script with the extracted values
    echo "Running with ws1=$ws1, ws2=$ws2, pmax1=$pmax1, pmax2=$pmax2, fout_name=$fout_name"
    ./run_model.jl "$ws1" "$ws2"  "$pmax1"  "$pmax2" "$fout_name" &  
    sleep 0.5   
done

echo "Alright, all jobs are submitted..."
#sleep for 10 min 
sleep 40m
echo "All runs are done!" 