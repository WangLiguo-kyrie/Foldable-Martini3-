**Tutorial:**

1. Polyply generate the topology without backbone dipole moment. 
`polyply gen_params` is used to generate the topology with protein sequence fasta input by `-seqf`(or protein sequence string by `-seq`). `-f` flag provides the force field. `polyply gen_coords` is used to generate a random IDR coordinates file by random walking with only IDR topology provided (`-p` flag). `-box` here defines the box of coordinates, `-o` defines the coordinates file name.
```
polyply gen_params -name p53TAD -f aminoacids-FoldableMartini3.ff -seqf p53.fasta -o p53TAD.itp
```
2. add_dipole.py add backbone dipole.
```
python add_dipole.py p53TAD
```

**Feature**

***Backbone dipole monment representation***

- Partial charge: 0.34e * 0.14nm
- More detailed LJ potential between UN-UP: based on helix formation penalty of amino acids, three types were classified.
  - UP1/UN1: A/R/L/K/M
  - UP2/UN2: E/Q/I/F/S/Y/W
  - UP3/UN3:N/D/C/H/T/V
  - Pro And Gly: Dipole charge and LJ protential were both removed for Pro (LJ for negative potetnial remained), only LJ potential was removed for Gly to reduce these two residue's ss formation propensity. 
- sigma=0.196,  epsilon: kJ/mol
  - UP1-UN1=3
  - UP1-UN2=2.5
  - UP1-UN3=2
  - UP2-UN1=2.5
  - UP2-UN2=2
  - UP2-UN3=1.5
  - UP3-UN1=2
  - UP3-UN2=1.5
  - UP3-UN3=1 
- dipole moment configuration:
  - orientation restrain occurs in the peptide bond planar, restrain UN in resid-i and UP in next resid-i+1 to be in opposite direction;
  - consider dipole shouldn't left/right shift just like sidechain, also add SBB/BBS angle.
 
***Backbone parameters***

- Multiple-minimas BBBB dihedral potential from ProMPT;
- BB-BB(-1)-BB(+1)-SC improper should be consistent in coil/Helix/sheet, so add this definition as in Martini-IDP to describe sidechain orientation.
- Pro specific BBBB dihedral for XXPX position to describe the helix disruption.
