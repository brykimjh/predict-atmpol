# predict-atmpol
Fortran code, which calculates the resulting gradients for molecular dynamics from the chaperone polarizabilities.

"Doubly Polarized QM/MM with Machine Learning Chaperone Polarizability," Kim, B.; Shao, Y.; Pu, J. J. Chem. Theory Comput. 2021, 17, 7682-7695 (doi:10.1021/acs.jctc.1c00567; PMID:34723536).  

Author: Bryant Kim  
Email: brykimjh@gmail.com  
Website: https://wikipugr.sitehost.iu.edu/wiki/index.php?title=The_Pu_Research_Group  

To compile test.f90 to user.exe and run user.exe (creates output.dat) run the following bash script:
```
$ bash 0run.sh
```

| File/Directory| Description   |
| ------------- | ------------- |
| ann.dat | artificial neural network which predicts atomic polarizabilities from geometry  |
| crd.txt | coordinates of qm/mm system  |
| force.txt | snapshot forces for crd.txt from a single step in molecular dynamics  |
| output.dat | output file containing updated forces from the chaperone polarizability corrections  |
| test.f90 | source code for dp-QM/MM correction  |
| user.exe | executeable for fortran code  |
