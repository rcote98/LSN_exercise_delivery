
echo "Launching instance of the simulation with the Verlet algorithm..."
echo

./../esercitazione4/simulate.sh

echo 
echo "DONE!"
echo 
echo "Copying the data files..."
echo

cp /../esercitazione4/sim_data/solid/output.gofr_ave.dat sim_data/solid/output.gofr-verlet.0
cp /../esercitazione4/sim_data/liquid/output.gofr_ave.dat sim_data/liquid/output.gofr-verlet.0
cp /../esercitazione4/sim_data/gas/output.gofr_ave.dat sim_data/gas/output.gofr-verlet.0

echo 
echo "DONE!"
echo 
echo "Please, continue your I-Python notebook activity."