# LSN Exercise Delivery

Repository created to deliver the exercises for the subject "Laboratorio di Simulazione Numerica" of the Physics Department in the University of Milan.

## Requirements

- g++ C++ compiler
- Python 3.* with numpy, scipy and maplotlib modules installed
- IPython Notebook

## Usage

In general, each folder contains the cpp code plus two IPython notebooks, one with the exercises named "LSN_Exercises_##" plus one with the solutions called "esercitazione#", where # is the exercise number.

To view the solutions, compile the cpp code and execute all the .exe files generated with

  cd [exercise session directory]
  make
  ./esercizio.\*.exe

The raw data will be stored in conveniently named .csv files. To view them, execute all the commands within the notebook in order (many variable names are reused for convenience).
