# CAPSens_design

This repository describes methods and provides examples for modeling receptor:peptide complexes and performing subsequent computational design and model refinement.

Requirements:

* [Python3](https://www.python.org/downloads/)
* [Rosetta](https://new.rosettacommons.org/demos/latest/tutorials/install_build/install_build)

Read the paper: [Jefferson RE et al, 2022](https://doi.org/10.1101/2022.03.30.486413)

# Receptor:Peptide Complex Modeling

Modeling of receptor:peptide complexes follows 3 overarching steps:

**1. Receptor scaffold hybridization**  
**1. Flexible peptide docking and diversification**  
**1. Loop modeling and structure relaxation**  

How to set up each step will depend on your specific target receptor:peptide complex, but guidelines are discussed below. The demo provides inputs for the CXCR4:CXCL12 receptor:peptide pair used for the generation of CAPSens designs. See `demo` subdirectories for specific instructions on how to setup and run simulations. The `scripts` folder has useful tools for running the protocol. The `CAPSen_models` contains final models generated by the protocol and molecular dynamics simulations.

## 1. Receptor scaffold hybridization

### Homologous template selection

Because active-state structures of membrane receptors is not always available, this method uses a homology modeling approach to generate possible scaffolds. We use [HHpred](https://toolkit.tuebingen.mpg.de/tools/hhpred) to identify homologous template structures. Homologous templates should share at least 20% sequence identity. Keep in mind that alignments may need to be manually adjusted to fix local mismatches, especially around intra-receptor disulfide bonds. The native sequence can be threaded onto the homologous template using Step 1 of the [RosettaCM protocol](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/RosettaCM). This uses the grishin alignment format and Rosetta's partial_thread binary.

### Hybridization

For hybridization, the structural elements to recombine should be selected with care. While a homologous template may be the best overall match to your target complex, they may be locally divergent segments (>10 residues) that are expected to make key contacts with the peptide ligand, where another template may have higher sequence identity, or where the native inactive structure may serve better to model the complex. The aim is to incorporate the maximal number of active-state features from the template while preventing significant _de novo_ reconstruction of the transmembrane core region due to poor sequence-structure alignment. It is good practice to use all combinations of structural elements to cover as much conformational space as possible, as the structure of the receptor scaffold can have significant effects in subsequent modeling steps.

In the case of CXCR4, the active-state template US28 (PDB ID: 4XT1) was used with 29% sequence identity. In the absence of other homologous active-state structures, the inactive-state structure of CXCR4 (PDB ID: 4RWS) was used to model local regions that differed significantly between both templates: ECL2 (residues 87-101) and the extracellular head of transmembrane helix (TM) 2 (residues 174-192). ECL2 of native CXCR4 forms a long beta-hairpin structure, which is shorter and not well-resolved in US28. The tip of TM2 differs due to a proline kink in the native sequence, which is not present in US28.

The compose script assembles hybrid scaffolds from donor and host pdbs. The example takes the extracellular head and ECL2 of CXCR4 (PDB ID: 4RWS) and places them into the homology template based on US28 (PDB ID: 4XT1).

```
cd 1_receptor-peptide_modeling/0_scaffold_hybridization/
../../../scripts/compose -d input/donor.pdb -f input/composition_file input/host.pdb >CXCR4.TM2-ECL2_4xt1_template.pdb
```

### Peptide threading

If the active-state template is complexed with a peptide, this can serve as a viable template for threading your peptide ligand of interest. Alternately, another complexed peptide template may serve better. Additionally, if a large enough pool of _apo_ peptide structures are available, you may find a consensus peptide conformation that could be used as well. This threading is an initial conformation that will be allowed all possible degrees of freedom in the following peptide docking step.

In the case of CXCR4, the US28 active-state structure was complexed with the CX3CL3 chemokine. CXCL12 was threaded onto the complexed peptide structure with K1 of CXCL12 aligned to the H2 position of CX3CL1 to match the partial positive charge of the imidazole ring. Q1 of CX3CL1 is cyclized to form pyroglutamate to produce a neutral N-terminus. Alternately, K1 of CXCL12 was aligned to H3 of CX3CL1, but docking from this initial position yielded models with weak interface energies and few contacts to key binding residues.

## 2. Flexible peptide docking and diversification

The peptide is translated and rotated across the binding pocket to saturate all possible binding modes. From these initial positions, the peptide ligand is docked with all degrees of flexibility and receptor repacking allowed at the interface. Optionally, experimentally informed constraints can be applied in this step fto favor putative interfacial contacts. The top-scoring decoys are filtered and diversified to select many different peptide poses by position, shape, and orientation. 200-220 unique peptide conformations are taken as inputs for the next step.

The pep_grid.py script generates inputs for flexible peptide docking by translating and rotating the peptide across the receptor binding pocket. You will probably want to adjust parameters for your particular system.

```
cd 1_receptor-peptide_modeling/1_peptide_docking/prep
../../../../scripts/pep_grid.py CXCR4-CXCL12_TM2.ECL2_hybrid.pdb
```

Prepack the inputs before running flexible peptide docking.

```
for n in {0..116}; do
    printf -v sn "%03d" $n
    $ROS/source/bin/FlexPepDocking.linuxgccrelease -in:file:spanfile ../input/CXCR4_SDF1.span -database $ROS/database/ -s $sn.pdb
done
```

## 3. Loop modeling and structure relaxation

Missing receptor loop residues are rebuilt _de novo_ around each diversified peptide pose and the receptor:peptide complexes are relaxed to simulate induced fit effects. Inter-TM constraints derived from sequence conservation are applied to restrain receptor structure. Any additional experimentally informed constraints used in peptide docking (step 2) can be applied again to preserve putative interfacial contacts during the full structure relax.

The 10 % top-scoring decoys are clustered by structural similarity. Representative models from clusters of sufficient size were analyzed for total score, buried surface area, interface score, peptide score, and key contacts with residues known to be important for receptor activation.

# Computational Design

Combinatorial mutagenesis using flexible backbone repacking.

If you have the means, it may be advantageous to test a small rational library of point mutants of receptor residues that have any heteroatom within 4 A of all selected peptide models.

# Model Refinement

The models that best reflect experimentally validated shifts in activity are re-docked in mutational contexts to explore receptor:peptide interactions that may not have been fully sampled in the modeling of the native WT complex. A single constraint was used to anchor a key electrostatic interaction in the depth of the binding pocket.
