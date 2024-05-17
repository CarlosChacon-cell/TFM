

# IMMUNOGENICITY 

In this notebook I will put everything I obtain from the immunogenicity study. 

As summary, immunogenicity is pretty pretty complicated, and probably having a global perspective is quite complex (Most likely impossible for me). However, I am going to try several approaches 

## B-CELL LINEAR EPITOPES

First, we start with the **antigenic_response.py** . This is a pretty simplistic method to compute the presence of linear epitopes of B-cell lymphocites. It is based on the propensity with which aminoacid are present in the surface and from there it computes the immunogenicity based on the frequency at which each residue is present in known epitopes.

Based on this very doubtful [webpage](http://imed.med.ucm.es/Tools/antigenic.pl), which is at the same time based on this old [paper](https://pubmed.ncbi.nlm.nih.gov/1702393/), I have measured the antigenicity of all the micropeptides.
I also made a bias, in which each aminoacid is ponderated by **1-aminocid_propensity** (negative for antigenic residues, positive for non-antigenic residues)

Then, I re-ran the micropeptide binders and compute if the antigenicity increased or decreased using the bias.

![alt text](figures/antigenic_propensity.png)

As can be seen in the plot, the bias decreases the antigenicity in almost all cases a mean value of 
-0.004765 while the positive bias increases the antigenic value a mean of 0.012092

Maybe multiplicating this bias by 2,3 or more could help

## MHC-I EPITOPES

Adaptative immunogenicity (the one we are interested in) can be elicited through two different sources, B-cells and T-cells. Inside T-cells we have CD4 (helpers) and CD8(killers). CD4 increase an immunogenic response when they are activated through the binding with HLA-II present in Antigen Presenting Cells (APC) These HLAs generally present peptides of exogenous proteins that the APCs have taken through endocytosis. On the other hand, C8 are activated through HLA-I, which bind proteins that have been tore down by the proteasome (so proteins present in the cytosol). Once  a CD8 is activated, it kills the cell presenting the antigen.

I use the [PRIME](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9811684/pdf/main.pdf#appsec1) code to compute the possible presence of peptides that can elicit the immunogenic response through binding to HLA-I. Following what they describe in the paper, we would consider as possible candidates those peptides whose %rank is smaller than 0.5. We would only use 9-mer peptides since most of the peptides has that lenght and using other lengths will complicate things (These are preliminary studies). Those peptides can be generated unsing **mhc-generator.py** 

The problem is that there are a lot of HLA alleles (PRIME gives you the possibility that the peptide binds any of the HLA alleles, detailing to which allele it will bind), so probably this computation can be used to select sequences that reduces the probability that the protein is recognized by the most frequent HLA (personalized medicine), but we cannot assure it will not elicit an immunogenic response.

I made a list with the most frequents alleles from type A, B and C in a European population, which is saved in /emdata/cchacon/inmuno/micropeptides/run_33_europe.txt

