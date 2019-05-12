#/usr/bin/bash

echo "Cleaning stuff up..."
echo " "

make clean
./clean.sh

echo " "
echo "Compiling..."
echo " "

make 

echo " "
echo "Starting Simulation Process..."

# SOLID SIMULATION

echo " "
echo "Starting SOLID simulation:"
echo " -- Equilibration phase -- "
echo " "

cp setups/input.solid input.dat
cp setups/config.fcc config.0

./Monte_Carlo_NVT.exe

mkdir sim_data/solid
mv output*.0 sim_data/solid/



echo " "
echo " -- Measurement Phase -- "
echo " "

cp config.final config.0

./Monte_Carlo_NVT.exe

mv output*.0 sim_data/solid/


echo " "
echo "DONE!"
echo " "

# LIQUID SIMULATION

echo " "
echo "Starting LIQUID simulation:"
echo " -- Equilibration phase -- "
echo " "


cp setups/input.liquid input.dat
cp setups/config.fcc config.0

./Monte_Carlo_NVT.exe

mkdir sim_data/liquid
mv output*.0 sim_data/liquid/


echo " "
echo "Yet again..."
echo " "

cp config.final config.0

./Monte_Carlo_NVT.exe

mv output*.0 sim_data/liquid/

echo " "
echo " -- Measurement Phase -- "
echo " "

cp config.final config.0

./Monte_Carlo_NVT.exe

mv output*.0 sim_data/liquid/

echo " "
echo "DONE!"
echo " "

# GAS SIMULATION

echo " "
echo "Starting GAS simulation:"
echo " "

cp setups/input.gas input.dat
cp setups/config.fcc config.0

./Monte_Carlo_NVT.exe

mkdir sim_data/gas
mv output*.0 sim_data/gas/

echo " "
echo "Yet again..."
echo " "

cp config.final config.0

./Monte_Carlo_NVT.exe

mv output*.0 sim_data/gas/

echo " "
echo " -- Measurement Phase -- "
echo " "

cp config.final config.0

./Monte_Carlo_NVT.exe

mv output*.0 sim_data/gas/



echo " "
echo "DONE!"
echo "Simulation Process Finished!"