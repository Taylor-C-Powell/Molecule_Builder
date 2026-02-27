# MolBuilder Examples

Self-contained scripts demonstrating the MolBuilder pipeline.

## SMILES to Manufacturing

Run the full process chemistry pipeline from a SMILES string:

```bash
# Default: acetic acid
python examples/smiles_to_manufacturing.py

# Ethanol
python examples/smiles_to_manufacturing.py CCO

# Benzene
python examples/smiles_to_manufacturing.py "c1ccccc1"

# Aspirin
python examples/smiles_to_manufacturing.py "CC(=O)Oc1ccccc1C(=O)O"
```

**Output includes:**

1. Molecule parsing and functional group detection
2. Retrosynthetic route planning
3. Step-by-step synthesis report
4. Cost breakdown at 1 kg scale
5. Safety and hazard assessment

## Requirements

```bash
pip install molbuilder
```
