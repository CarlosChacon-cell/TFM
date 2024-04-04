

Este cuaderno esta destinado a ordenar un poco las cosas de cara a mi TFM

TFM
pMPNN improvement of TFEB
interaction scoring: write down and think
interaction scoring: compare with others (silent?)
Describe what you have done with TFEB and discuss
RFD
gather all data from previous runs
run interaction scoring on all and see if trends are similar
new designs?
improve new designs iteratively with pMPNN?
MD
Run with proteins
Run with small molecule hits
Implement short protocol to evaluate binding
Small Molecule
Implement RF-AllAtom and evaluate with ETP dataset
If better, redo, if worse, go back to think about clustering, etc.

Estas son las 4 patas del proyecto ahora mismo

## PMPNN

Esta sección ya esta terminada, no creo que haya más código que escribir.

Se ejecuta mediante:

``` bash 
/home/cchacon/Carlos_scripts/pmpnn_cycle_not_update.sh --input example.pdb --indices 2 3 4 5 27-50 --max 10 --peptide A --chains B
```

Donde input es el pdb con la cadena a modificar y el target, indices son los aminoacidos a fijar, max es el numero de secuencias que se van a crear, peptide es la cadena del binder y target la cadena del target.