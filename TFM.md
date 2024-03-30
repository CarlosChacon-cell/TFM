

## INTRODUCTION

This is intended to be something like a lab notebook in which I'm going to save all the results relevant for my TFM.

My TFM has 3 different legs (**PMPNN**, **RFD**,**RFAA** and **GROMACS**), so a different chapter will be dedicated to each. All of this computational methods will be employed with an example problem: RagAC and TFEB. So I'll start explaining this problem

## RagAC and TFEB.

RagA and RagC are two monomers that take place in the mTORC1-TFEB-Rag-Ragulator megacomplex (PDB: 7UXC).

TFEB is a transcription factor that act as master regulator of lysosomal biogenesis and autophagy. It is a bHLH-Leu zipper transcription factor whose activation increases the autophagic and lysosomal clearance capacity. Under normal circumstances, TFEB is regulated by cellular nutrient status through the phosphorylation of serines (Ser122/142/211) by mTORC1 under nutrient-replete conditions, allowing TFEB cytosolic retention and inactivation. To perform this inactivation, mTORC1 and TFEB interact through the RagGTPases, specially RagC GDP

TFEB phosphorylation is strictly dependent on the activation of FLCN (RagC GTPase activating protein), which remains inactive under aminoacid starvation.

The complex RagA-RagC-TFEB-Raptor was obtained by coexpression and its structure solve:

![alt text](/figures/complex_structure.png "Complex structure")

They also perform different mutation assays that help them gaina different insights into TFEB inactivation.

The **S211A** mutation stabilizes the TFEB association with the RagGTPases in cells.

**TFEB residues 2-105 are both necessary and sufficient** to form a stable complex with the RagGTPases: The first 40 residues form an alpha helix that occupies the cleft between the two G domains of the non-canonical (nc) dimer of RagA GTP and RagC GDP. Residues 2-18 of TFEB form a 570 $$1\overset{\circ}{A}$$ &#8491  interface