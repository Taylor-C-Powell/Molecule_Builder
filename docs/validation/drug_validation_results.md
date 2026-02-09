# MolBuilder Drug Validation Report

**Generated**: 2026-02-09 06:12:05
**MolBuilder Version**: 1.1.1
**Drugs Tested**: 50

## Summary

| Metric | Value |
|--------|-------|
| Total drugs | 50 |
| SMILES parsed | 50/50 (100%) |
| Retrosynthesis tree built | 50/50 (100%) |
| Disconnections found | 50/50 (100%) |
| Functional groups detected | 50/50 (100%) |
| Route coverage (non-macrocycle) | 100% |
| Errors | 0 |
| Total runtime | 79.7s |
| Average per drug | 1.59s |

## Per-Drug Results

| Drug | Parsed | Tree | FGs | Disconnections | Time (s) | Notes |
|------|--------|------|-----|----------------|----------|-------|
| aspirin | Y | Y | 4 | 5 | 0.44 | - |
| ibuprofen | Y | Y | 3 | 5 | 0.46 | - |
| acetaminophen | Y | Y | 3 | 5 | 0.35 | - |
| caffeine | Y | Y | 6 | 5 | 0.50 | - |
| naproxen | Y | Y | 5 | 5 | 0.77 | - |
| atorvastatin | Y | Y | 11 | 5 | 3.94 | - |
| metformin | Y | Y | 5 | 5 | 0.10 | - |
| omeprazole | Y | Y | 9 | 5 | 0.51 | - |
| losartan | Y | Y | 8 | 5 | 2.11 | - |
| lisinopril | Y | Y | 8 | 5 | 2.31 | - |
| amoxicillin | Y | Y | 7 | 5 | 2.08 | - |
| ciprofloxacin | Y | Y | 9 | 5 | 2.21 | - |
| sildenafil | Y | Y | 8 | 5 | 2.74 | - |
| sertraline | Y | Y | 5 | 5 | 0.50 | - |
| warfarin | Y | Y | 6 | 5 | 2.57 | - |
| fluoxetine | Y | Y | 7 | 5 | 0.92 | - |
| diazepam | Y | Y | 5 | 5 | 0.56 | - |
| gabapentin | Y | Y | 3 | 5 | 0.67 | - |
| metoprolol | Y | Y | 5 | 5 | 1.57 | - |
| amlodipine | Y | Y | 9 | 5 | 2.80 | - |
| furosemide | Y | Y | 8 | 5 | 0.95 | - |
| cetirizine | Y | Y | 8 | 5 | 2.56 | - |
| loratadine | Y | Y | 7 | 5 | 2.17 | - |
| captopril | Y | Y | 4 | 5 | 0.37 | - |
| enalapril | Y | Y | 6 | 5 | 1.54 | - |
| valsartan | Y | Y | 7 | 5 | 2.16 | - |
| methotrexate | Y | Y | 15 | 5 | 3.12 | - |
| tamoxifen | Y | Y | 6 | 5 | 1.02 | - |
| prednisone | Y | Y | 7 | 5 | 2.50 | - |
| dexamethasone | Y | Y | 8 | 5 | 2.76 | - |
| citalopram | Y | Y | 6 | 5 | 0.44 | - |
| venlafaxine | Y | Y | 4 | 5 | 1.30 | - |
| bupropion | Y | Y | 4 | 5 | 0.73 | - |
| lamotrigine | Y | Y | 7 | 5 | 0.75 | - |
| montelukast | Y | Y | 10 | 5 | 3.96 | - |
| pantoprazole | Y | Y | 12 | 5 | 0.43 | - |
| duloxetine | Y | Y | 5 | 5 | 0.69 | - |
| aripiprazole | Y | Y | 8 | 5 | 1.35 | - |
| rosuvastatin | Y | Y | 11 | 5 | 2.01 | - |
| esomeprazole | Y | Y | 9 | 5 | 0.49 | - |
| clopidogrel | Y | Y | 5 | 5 | 0.49 | - |
| levothyroxine | Y | Y | 11 | 5 | 0.98 | - |
| allopurinol | Y | Y | 4 | 5 | 0.14 | - |
| ranitidine | Y | Y | 7 | 5 | 0.90 | - |
| simvastatin | Y | Y | 6 | 5 | 3.97 | - |
| propranolol | Y | Y | 5 | 5 | 1.25 | - |
| methylphenidate | Y | Y | 4 | 5 | 0.87 | - |
| clonazepam | Y | Y | 6 | 5 | 0.57 | - |
| azithromycin | Y | Y | 13 | 5 | 7.41 | macrocycle |
| doxycycline | Y | Y | 10 | 5 | 3.71 | macrocycle |

## Methodology

Each drug's canonical SMILES was parsed by MolBuilder's SMILES parser, then fed into the retrosynthetic analysis engine with `max_depth=5` and `beam_width=5`. The engine identifies functional groups, proposes disconnections using the reaction knowledge base (181 templates), and builds a retrosynthetic tree toward purchasable starting materials (250+ entries).

Complex macrocycles (azithromycin, doxycycline) are expected to be challenging and are evaluated on whether the analysis completes without errors rather than on route quality.

## Compliance Notes

This validation report supports IQ/OQ documentation for 21 CFR Part 11 compliance. The test suite is fully automated and reproducible via:

```bash
python -m pytest tests/test_drug_validation.py -v
python tests/validation_report.py
```
