"""Toolbar: element palette, bond tools, mode buttons."""

import tkinter as tk
from tkinter import ttk


class MolToolbar(ttk.Frame):
    """Toolbar with element selector, bond tools, and action buttons."""

    COMMON_ELEMENTS = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"]
    BOND_TYPES = [("Single", 1), ("Double", 2), ("Triple", 3)]

    def __init__(self, parent, on_element_selected=None, on_bond_selected=None,
                 on_action=None):
        super().__init__(parent)
        self._on_element = on_element_selected
        self._on_bond = on_bond_selected
        self._on_action = on_action
        self.selected_element = tk.StringVar(value="C")
        self.selected_bond = tk.IntVar(value=1)
        self._build()

    def _build(self):
        # Element palette
        ef = ttk.LabelFrame(self, text="Element")
        ef.pack(side="left", padx=5, pady=2)
        for sym in self.COMMON_ELEMENTS:
            b = ttk.Radiobutton(ef, text=sym, value=sym,
                                variable=self.selected_element,
                                command=self._elem_changed)
            b.pack(side="left", padx=1)

        # Custom element entry
        self._custom = ttk.Entry(ef, width=4)
        self._custom.pack(side="left", padx=2)
        self._custom.bind("<Return>", self._custom_elem)

        # Bond type
        bf = ttk.LabelFrame(self, text="Bond")
        bf.pack(side="left", padx=5, pady=2)
        for label, val in self.BOND_TYPES:
            b = ttk.Radiobutton(bf, text=label, value=val,
                                variable=self.selected_bond,
                                command=self._bond_changed)
            b.pack(side="left", padx=1)

        # Action buttons
        af = ttk.LabelFrame(self, text="Actions")
        af.pack(side="left", padx=5, pady=2)
        for label in ["Add Atom", "Add Bond", "Delete", "Add H", "Clear"]:
            b = ttk.Button(af, text=label,
                           command=lambda l=label: self._action(l))
            b.pack(side="left", padx=2)

    def _elem_changed(self):
        if self._on_element:
            self._on_element(self.selected_element.get())

    def _custom_elem(self, event):
        sym = self._custom.get().strip().capitalize()
        if sym:
            self.selected_element.set(sym)
            self._elem_changed()

    def _bond_changed(self):
        if self._on_bond:
            self._on_bond(self.selected_bond.get())

    def _action(self, action):
        if self._on_action:
            self._on_action(action)
