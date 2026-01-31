"""Sidebar: molecule properties, analysis options."""

import tkinter as tk
from tkinter import ttk


class MolSidebar(ttk.Frame):
    """Side panel showing molecule properties and analysis controls."""

    def __init__(self, parent, on_analyze=None):
        super().__init__(parent, width=280)
        self.pack_propagate(False)
        self._on_analyze = on_analyze
        self._build()

    def _build(self):
        # Info panel
        info = ttk.LabelFrame(self, text="Molecule Info")
        info.pack(fill="x", padx=5, pady=5)

        self.name_var = tk.StringVar(value="(none)")
        self.formula_var = tk.StringVar(value="")
        self.atoms_var = tk.StringVar(value="0 atoms")
        self.bonds_var = tk.StringVar(value="0 bonds")
        self.mw_var = tk.StringVar(value="MW: --")

        for var in [self.name_var, self.formula_var, self.atoms_var,
                    self.bonds_var, self.mw_var]:
            ttk.Label(info, textvariable=var).pack(anchor="w", padx=5)

        # Analysis buttons
        af = ttk.LabelFrame(self, text="Analysis")
        af.pack(fill="x", padx=5, pady=5)
        analyses = [
            "Functional Groups",
            "Bond Analysis",
            "Stereochemistry",
            "SMILES",
            "Retrosynthesis",
            "Process Engineering",
        ]
        for name in analyses:
            b = ttk.Button(af, text=name,
                           command=lambda n=name: self._analyze(n))
            b.pack(fill="x", padx=5, pady=1)

        # Export buttons
        ef = ttk.LabelFrame(self, text="Export")
        ef.pack(fill="x", padx=5, pady=5)
        for fmt in ["XYZ", "MOL/SDF", "PDB", "JSON", "SMILES"]:
            b = ttk.Button(ef, text=f"Save as {fmt}",
                           command=lambda f=fmt: self._analyze(f"export_{f}"))
            b.pack(fill="x", padx=5, pady=1)

        # Results text area
        rf = ttk.LabelFrame(self, text="Results")
        rf.pack(fill="both", expand=True, padx=5, pady=5)
        self.results_text = tk.Text(rf, wrap="word", bg="#111122",
                                    fg="#ccddee", font=("Consolas", 9),
                                    height=12)
        scroll = ttk.Scrollbar(rf, orient="vertical",
                               command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        self.results_text.pack(fill="both", expand=True)

    def update_info(self, mol):
        """Update displayed molecule info."""
        if mol is None:
            self.name_var.set("(none)")
            self.formula_var.set("")
            self.atoms_var.set("0 atoms")
            self.bonds_var.set("0 bonds")
            self.mw_var.set("MW: --")
            return

        self.name_var.set(mol.name or "(unnamed)")

        # Compute formula
        from collections import Counter
        counts = Counter(a.symbol for a in mol.atoms)
        formula = ""
        for sym in ["C", "H", "N", "O", "S", "P"]:
            if sym in counts:
                formula += sym + (str(counts[sym]) if counts[sym] > 1 else "")
                del counts[sym]
        for sym in sorted(counts):
            formula += sym + (str(counts[sym]) if counts[sym] > 1 else "")
        self.formula_var.set(formula)

        self.atoms_var.set(f"{len(mol.atoms)} atoms")
        self.bonds_var.set(f"{len(mol.bonds)} bonds")

        from molbuilder.core.elements import atomic_weight
        mw = sum(atomic_weight(a.symbol) for a in mol.atoms)
        self.mw_var.set(f"MW: {mw:.2f} g/mol")

    def show_results(self, text: str):
        """Display analysis results."""
        self.results_text.delete("1.0", "end")
        self.results_text.insert("1.0", text)

    def _analyze(self, name):
        if self._on_analyze:
            self._on_analyze(name)
