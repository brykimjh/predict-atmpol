# predict-atmpol
Fortran code, which calculates the resulting gradients for molecular dynamics from the chaperone polarizabilities.

$ bash 0run.sh to compiles test.f90 to user.exe and run user.exe (creates output.dat)

| File/Directory| Description   |
| ------------- | ------------- |
| ann.dat | artificial neural network which predicts atomic polarizabilities from geometry  |
| crd.txt | coordinates of qm/mm system  |
| force.txt | snapshot forces for crd.txt from a single step in molecular dynamics  |
| output.dat | output file containing updated forces from the chaperone polarizability corrections  |
| test.f90 | source code for dp-QM/MM correction  |
| user.exe | executeable for fortran code  |
