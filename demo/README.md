# Running the Demo

## Receptor-Peptide Modeling

### Hybridization of receptor scaffold

Once you have set up appropriate homology templates and structures to use for scaffold hybridization, align them to each other and decide which local regions to hybridize. Once prepared you may use the compose script to take structural elements defined in a composition file from a donor pdb and stitch them into a host pdb. The example takes the extracellular head and ECL2 of CXCR4 (PDB ID: 4RWS) and places them into the homology template based on US28 (PDB ID: 4XT1).

```
cd 1_receptor-peptide_modeling/0_scaffold_hybridization/
../../../scripts/compose -d input/donor.pdb -f input/composition_file input/host.pdb >CXCR4.TM2-ECL2_4xt1_template.pdb
```

### Flexible peptide docking

Depending on active-state structures with complexed peptide ligands available for homologous templates, you can try threading your target on several peptide structures. In the example, CXCL12 N-terminal sequence was threaded onto the CX3CL1