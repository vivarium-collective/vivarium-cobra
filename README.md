# Vivarium-COBRA

A vivarium wrapper process for the COBRApy library.
`COBRA_FBA` can be loaded with BiGG models, and has been primarily used with BiGG model iAF1260b. 
The project includes `CobraComposite`, which adds auxiliary processes `MassDeriver` and `VolumeDeriver` 
to make a type of dynamic FBA. Other processes can constrain individual fluxes in `COBRA_FBA` through 
the flux_bounds port, allowing for a type of integrative FBA.

## Installation

```
$ pip install -r requirements.txt
```
