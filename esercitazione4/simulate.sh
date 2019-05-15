#/usr/bin/bash

echo "Cleaning stuff up..."
echo

make clean
./clean.sh

echo 
echo "Compiling..."
echo 

make 

echo 
echo "Starting Simulation Process..."
echo

mkdir -p sim_data

echo "------------------------------------------------------"
# SOLID SIMULATION--------------------------------------

echo 
echo "Starting SOLID simulation:"
echo
echo " ### Equilibration phase ### "
echo 

cp setups/input.solid.restart input.dat
cp setups/config.fcc config.0

./esercizio4.exe

mkdir -p sim_data/solid
mv output.*.dat sim_data/solid/


echo 
echo " ### Measurement Phase ### "
echo 

cp setups/input.solid input.dat
cp config.final config.0
cp config.final.prev config.0.prev

./esercizio4.exe

mv output.*.dat sim_data/solid/
mv config.final sim_data/solid/
mv config.final.prev sim_data/solid/

echo 
echo "DONE!"
echo 

echo "------------------------------------------------------"
# LIQUID SIMULATION ------------------------------------

echo 
echo "Starting LIQUID simulation:"
echo
echo " ### Equilibration phase ### "
echo 

cp setups/input.liquid.restart input.dat
cp setups/config.fcc config.0



./esercizio4.exe

mkdir -p sim_data/liquid
mv output.*.dat sim_data/liquid/


echo 
echo "Yet again..."
echo 

cp setups/input.liquid input.dat

cp config.final config.0
cp config.final.prev config.0.prev

./esercizio4.exe

mv output.*.dat sim_data/liquid/

echo 
echo " ### Measurement Phase ### "
echo 

cp config.final config.0
cp config.final.prev config.0.prev

./esercizio4.exe

mv output.*.dat sim_data/liquid/
mv config.final sim_data/solid/
mv config.final.prev sim_data/solid/

echo
echo "DONE!"
echo

echo "------------------------------------------------------"
# GAS SIMULATION ---------------------------------------

echo
echo "Starting GAS simulation:"
echo
echo " ### Equilibration phase ### "
echo 

cp setups/input.gas.restart input.dat
cp setups/config.fcc config.0

./esercizio4.exe

mkdir -p sim_data/gas
mv output.*.dat sim_data/gas/

echo
echo "Yet again..."
echo

cp setups/input.gas.restart input.dat

cp config.final config.0
cp config.final.prev config.0.prev

./esercizio4.exe

mv output.*.dat sim_data/gas/

echo 
echo " ### Measurement Phase ### "
echo 

cp config.final config.0
cp config.final.prev config.0.prev

./esercizio4.exe

mv output.*.dat sim_data/gas/
mv config.final sim_data/solid/
mv config.final.prev sim_data/solid/

echo
echo "DONE!"

echo "------------------------------------------------------"

echo
echo "Simulation Process Finished!"
echo
echo "Restoring default simulation parameters..."

cp setups/input.def input.dat
cp setups/config.fcc config.0
rm config.0.prev

echo "DONE!"

echo
echo "SIMULATION PROCESS COMPLETE!"
echo "Now you can visualize the results on the I-Python Notebook."
