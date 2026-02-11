"""Dialog windows: SMILES input, file open/save, export settings."""

import tkinter as tk
from tkinter import ttk, filedialog


class SmilesDialog:
    """Dialog for entering a SMILES string."""

    def __init__(self, parent):
        self.result = None
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Enter SMILES")
        self.dialog.geometry("400x150")
        self.dialog.transient(parent)
        self.dialog.grab_set()

        ttk.Label(self.dialog, text="SMILES string:").pack(padx=10, pady=(10, 0), anchor="w")
        self.entry = ttk.Entry(self.dialog, width=50)
        self.entry.pack(padx=10, pady=5, fill="x")
        self.entry.focus_set()
        self.entry.bind("<Return>", lambda e: self._ok())

        ttk.Label(self.dialog, text="Name (optional):").pack(padx=10, anchor="w")
        self.name_entry = ttk.Entry(self.dialog, width=50)
        self.name_entry.pack(padx=10, pady=5, fill="x")

        bf = ttk.Frame(self.dialog)
        bf.pack(pady=10)
        ttk.Button(bf, text="OK", command=self._ok).pack(side="left", padx=5)
        ttk.Button(bf, text="Cancel", command=self._cancel).pack(side="left", padx=5)

        parent.wait_window(self.dialog)

    def _ok(self):
        smi = self.entry.get().strip()
        name = self.name_entry.get().strip()
        if smi:
            self.result = (smi, name)
        self.dialog.destroy()

    def _cancel(self):
        self.dialog.destroy()


class FormulaDialog:
    """Dialog for entering a molecular formula (for simple VSEPR molecules)."""

    def __init__(self, parent):
        self.result = None
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Enter Formula")
        self.dialog.geometry("350x120")
        self.dialog.transient(parent)
        self.dialog.grab_set()

        ttk.Label(self.dialog, text="Molecular formula (e.g. H2O, CH4):").pack(padx=10, pady=(10, 0), anchor="w")
        self.entry = ttk.Entry(self.dialog, width=30)
        self.entry.pack(padx=10, pady=5)
        self.entry.focus_set()
        self.entry.bind("<Return>", lambda e: self._ok())

        ttk.Label(self.dialog, text="Charge:").pack(padx=10, anchor="w")
        self.charge = ttk.Entry(self.dialog, width=5)
        self.charge.insert(0, "0")
        self.charge.pack(padx=10, pady=2, anchor="w")

        bf = ttk.Frame(self.dialog)
        bf.pack(pady=5)
        ttk.Button(bf, text="OK", command=self._ok).pack(side="left", padx=5)
        ttk.Button(bf, text="Cancel", command=self._cancel).pack(side="left", padx=5)

        parent.wait_window(self.dialog)

    def _ok(self):
        formula = self.entry.get().strip()
        try:
            charge = int(self.charge.get().strip())
        except ValueError:
            charge = 0
        if formula:
            self.result = (formula, charge)
        self.dialog.destroy()

    def _cancel(self):
        self.dialog.destroy()


def file_open_dialog(parent, filetypes=None):
    """Show file open dialog, return filepath or None."""
    if filetypes is None:
        filetypes = [
            ("All supported", "*.xyz *.mol *.sdf *.pdb *.json *.smi"),
            ("XYZ files", "*.xyz"),
            ("MOL files", "*.mol"),
            ("SDF files", "*.sdf"),
            ("PDB files", "*.pdb"),
            ("JSON files", "*.json"),
            ("SMILES files", "*.smi"),
            ("All files", "*.*"),
        ]
    return filedialog.askopenfilename(parent=parent, filetypes=filetypes)


def file_save_dialog(parent, default_ext=".xyz", filetypes=None):
    """Show file save dialog, return filepath or None."""
    if filetypes is None:
        filetypes = [
            ("XYZ files", "*.xyz"),
            ("MOL files", "*.mol"),
            ("SDF files", "*.sdf"),
            ("PDB files", "*.pdb"),
            ("JSON files", "*.json"),
            ("SMILES files", "*.smi"),
        ]
    return filedialog.asksaveasfilename(parent=parent, defaultextension=default_ext,
                                        filetypes=filetypes)
