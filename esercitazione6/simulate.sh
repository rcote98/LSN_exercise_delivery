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

mkdir -p "sim_data/"

REP_NUMBER=(50)

for (( c=0; c<=$REP_NUMBER; c++ ))
do  
    TEMP=`echo "scale=3; 0.5 + 1.5*$c/$REP_NUMBER" | bc`

    mkdir -p "sim_data/TEMP$c"

    echo "#### Simulation Round #$c of $REP_NUMBER #####"
    echo

    # ------------------------------------------------------

    echo "Starting Metropolis Simulations:"
    echo

    mkdir -p "sim_data/TEMP$c/metro"

    echo "WITHOUT magnetization:"
    echo "Equilibration..."
    echo

    cp setups/input.metro.eq input.dat

    ./esercizio6.exe $TEMP
    
    echo
    echo "Measurement..."
    echo

    cp setups/input.metro input.dat

    ./esercizio6.exe $TEMP

    mv output.*.0 "sim_data/TEMP$c/metro"
    mv conf.final "sim_data/TEMP$c/metro"

    echo "WITH magnetization:"
    echo "Equilibration..."
    echo

    cp setups/input.metro.eq input.dat

    ./esercizio6.exe $TEMP 0.02
    
    echo
    echo "Measurement..."
    echo

    cp setups/input.metro input.dat

    ./esercizio6.exe $TEMP 0.02

    mv output.mag.0 "sim_data/TEMP$c/metro"
    mv conf.final "sim_data/TEMP$c/metro/conf.final.mag"
    rm output.*.0

    echo


    # ------------------------------------------------------

    echo "Starting Gibbs Simulations:"
    echo

    mkdir -p "sim_data/TEMP$c/gibbs"

    echo "WITHOUT magnetization:"
    echo "Equilibration..."
    echo

    cp setups/input.gibbs.eq input.dat

    ./esercizio6.exe $TEMP
    
    echo
    echo "Measurement..."
    echo

    cp setups/input.gibbs input.dat

    ./esercizio6.exe $TEMP

    mv output.*.0 "sim_data/TEMP$c/gibbs"
    mv config.final "sim_data/TEMP$c/gibbs"

    echo "WITH magnetization:"
    echo "Equilibration..."
    echo

    cp setups/input.gibbs.eq input.dat

    ./esercizio6.exe $TEMP 0.02
    
    echo
    echo "Measurement..."
    echo

    cp setups/input.gibbs input.dat

    ./esercizio6.exe $TEMP 0.02

    mv output.mag.0 "sim_data/TEMP$c/gibbs"
    mv config.final   "sim_data/TEMP$c/gibbs/conf.final.mag"
    rm output.*.0

    echo 
    echo "DONE!"
    echo

done

echo "Restoring default simulation parameters..."

cp setups/input.def input.dat

echo "DONE!"

echo
echo "SIMULATION PROCESS COMPLETE!"
echo "Now you can visualize the results on the I-Python Notebook."
