

## INTRODUCTION

This is intended to be something like a lab notebook in which I'm going to save all the results relevant for my TFM.

My TFM has 3 different legs (**PMPNN**, **RFD**,**RFAA** and **GROMACS**), so a different chapter will be dedicated to each. All of this computational methods will be employed with an example problem: RagAC and TFEB. So I'll start explaining this problem

## RagAC and TFEB.

RagA and RagC are two monomers that take place in the mTORC1-TFEB-Rag-Ragulator megacomplex (PDB: 7UXC).

TFEB is a transcription factor that act as master regulator of lysosomal biogenesis and autophagy. It is a bHLH-Leu zipper transcription factor whose activation increases the autophagic and lysosomal clearance capacity. Under normal circumstances, TFEB is regulated by cellular nutrient status through the phosphorylation of serines (Ser122/142/211) by mTORC1 under nutrient-replete conditions, allowing TFEB cytosolic retention and inactivation. To perform this inactivation, mTORC1 and TFEB interact through the RagGTPases, specially RagC GDP

TFEB phosphorylation is strictly dependent on the activation of FLCN (RagC GTPase activating protein), which remains inactive under aminoacid starvation.

The complex RagA-RagC-TFEB-Raptor was obtained by coexpression and its structure solve:

![alt text](/figures/complex_structure.png "Complex structure")

### STRUCTURE DETAILS

The **S211A** mutation stabilizes the TFEB association with the RagGTPases in cells.

**TFEB residues 2-105 are both necessary and sufficient** to form a stable complex with the RagGTPases: **The first 40 residues form an alpha helix that occupies the cleft** between the two G domains of the non-canonical (nc) dimer of RagA GTP and RagC GDP. Residues 2-18 of TFEB form a 570 A interface with the RagC GDP G domain.

TFEB is clamped in the cleft by **hydrogen bonds between RagC D294 and the backbone of the TFEB R4 and I5**, and **salt bridges between RagC D290 and TFEB R4 and R8**. More hydrophobic contacts between TFEB L7/11 and RagC V80, I220 and between TFEB Q15 sidechain and RagC Y221.

TFEB interacts with RagA through a stacking interaction between the carbon chain of R13 and RagA W165 and hydrophobic contacts between TFEB I5, M9 nad RagA I234. Moreover, there are also another interaction between TFEB 76-83 and RagA GTP alpha 4 (TFEB L76 with a hydrophobic pocket and a salt bridge between TFEB H84 and RagA E111). M

**TFEB mutants (I5D, L7D/R8D and M9D) were not phosphorylated** as well as **RagA mutants (H104D, Q107D,E111R)**, and thus is localized in the nucleus.



RagC clamp is required for TFEB phosphorylation,and when we mutate D294R we avoid the phosphorylation of TFEB, which is then allocated in the nucleus.












05/04
Aqui voy a apuntar solo unas cosillas que si no se me olvidan:
Ya tenemos el pae_interaction local funcionando y un pipeline para correrlo sobre los hits, y voy a modificar el submit_af2_interfaces para que coja el cutre
Le he mandado a Rafa las mutantes de PolG a ver si que hay que hacer mas
En ashoka parece que el af2 tarda approx 100 segundos menos que en la maquina local
Working in ashoka and ready to go
