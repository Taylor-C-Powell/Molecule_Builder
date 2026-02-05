"""Main tkinter application for the MolBuilder 3D GUI."""

import tkinter as tk
from tkinter import ttk, messagebox

from molbuilder.gui.canvas3d import MolCanvas3D
from molbuilder.gui.toolbar import MolToolbar
from molbuilder.gui.sidebar import MolSidebar
from molbuilder.gui.dialogs import SmilesDialog, FormulaDialog, file_open_dialog, file_save_dialog
from molbuilder.gui.event_handler import MolEventHandler
from molbuilder.molecule.graph import Molecule


class MolBuilderApp:
    """Main MolBuilder GUI application."""

    def __init__(self):
        self.root = tk.Tk()
        self.root.title("MolBuilder - Molecular Engineering Tool")
        self.root.geometry("1200x800")
        self.root.configure(bg="#1a1a2e")

        self.handler = MolEventHandler()
        self._build_menu()
        self._build_ui()
        self._refresh()

    def _build_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="New", command=self._new)
        file_menu.add_command(label="Open...", command=self._open)
        file_menu.add_command(label="Save As...", command=self._save)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)

        # Build menu
        build_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Build", menu=build_menu)
        build_menu.add_command(label="From SMILES...", command=self._from_smiles)
        build_menu.add_command(label="From Formula (VSEPR)...", command=self._from_formula)
        build_menu.add_separator()
        build_menu.add_command(label="Preset: Ethanol", command=lambda: self._preset("CCO", "ethanol"))
        build_menu.add_command(label="Preset: Aspirin", command=lambda: self._preset("CC(=O)Oc1ccccc1C(=O)O", "aspirin"))
        build_menu.add_command(label="Preset: Caffeine", command=lambda: self._preset("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "caffeine"))
        build_menu.add_command(label="Preset: Benzene", command=lambda: self._preset("c1ccccc1", "benzene"))
        build_menu.add_command(label="Preset: Glucose", command=lambda: self._preset("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "glucose"))

        # Analysis menu
        analysis_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)
        analysis_menu.add_command(label="Detect Functional Groups", command=lambda: self._run_analysis("Functional Groups"))
        analysis_menu.add_command(label="Bond Analysis", command=lambda: self._run_analysis("Bond Analysis"))
        analysis_menu.add_command(label="Generate SMILES", command=lambda: self._run_analysis("SMILES"))

        # Simulate menu
        simulate_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Simulate", menu=simulate_menu)
        simulate_menu.add_command(
            label="MD Vibration...", command=self._sim_md_vibration)
        simulate_menu.add_command(
            label="Bond Formation...", command=self._sim_bond_formation)
        simulate_menu.add_command(
            label="SN2 Mechanism...", command=self._sim_sn2)
        simulate_menu.add_separator()
        simulate_menu.add_command(
            label="Export Animation...", command=self._sim_export)

    def _build_ui(self):
        # Toolbar at top
        self.toolbar = MolToolbar(
            self.root,
            on_element_selected=self._on_element,
            on_bond_selected=self._on_bond,
            on_action=self._on_action,
        )
        self.toolbar.pack(fill="x")

        # Main area: canvas + sidebar
        main = ttk.PanedWindow(self.root, orient="horizontal")
        main.pack(fill="both", expand=True)

        canvas_frame = ttk.Frame(main)
        self.canvas = MolCanvas3D(canvas_frame)
        main.add(canvas_frame, weight=3)

        self.sidebar = MolSidebar(main, on_analyze=self._run_analysis)
        main.add(self.sidebar, weight=1)

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status = ttk.Label(self.root, textvariable=self.status_var, relief="sunken")
        status.pack(fill="x", side="bottom")

    def _refresh(self):
        """Re-render molecule and update sidebar."""
        self.canvas.render(self.handler.mol)
        self.sidebar.update_info(self.handler.mol)

    def _on_element(self, sym):
        self.handler.current_element = sym
        self.status_var.set(f"Element: {sym}")

    def _on_bond(self, order):
        self.handler.current_bond_order = order
        labels = {1: "Single", 2: "Double", 3: "Triple"}
        self.status_var.set(f"Bond: {labels.get(order, '?')}")

    def _on_action(self, action):
        if action == "Add Atom":
            self.handler.add_atom_free()
        elif action == "Add Bond":
            if not self.handler.add_bond_between_selected():
                self.status_var.set("Select 2 atoms to bond")
                return
        elif action == "Delete":
            self.handler.delete_selected()
        elif action == "Add H":
            self.handler.add_hydrogens()
        elif action == "Clear":
            self.handler.clear()
        self._refresh()

    def _new(self):
        self.handler.clear()
        self.handler.mol.name = "untitled"
        self._refresh()

    def _open(self):
        path = file_open_dialog(self.root)
        if not path:
            return
        try:
            ext = path.rsplit(".", 1)[-1].lower()
            if ext == "xyz":
                from molbuilder.io.xyz import read_xyz
                self.handler.mol = read_xyz(path)
            elif ext in ("mol", "sdf"):
                from molbuilder.io.mol_sdf import read_mol
                self.handler.mol = read_mol(path)
            elif ext == "pdb":
                from molbuilder.io.pdb import read_pdb
                self.handler.mol = read_pdb(path)
            elif ext == "json":
                from molbuilder.io.json_io import read_json
                self.handler.mol = read_json(path)
            else:
                messagebox.showerror("Error", f"Unsupported format: .{ext}")
                return
            self.status_var.set(f"Opened: {path}")
            self._refresh()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _save(self):
        path = file_save_dialog(self.root)
        if not path:
            return
        try:
            ext = path.rsplit(".", 1)[-1].lower()
            if ext == "xyz":
                from molbuilder.io.xyz import write_xyz
                write_xyz(self.handler.mol, path)
            elif ext == "mol":
                from molbuilder.io.mol_sdf import write_mol
                write_mol(self.handler.mol, path)
            elif ext == "pdb":
                from molbuilder.io.pdb import write_pdb
                write_pdb(self.handler.mol, path)
            elif ext == "json":
                from molbuilder.io.json_io import write_json
                write_json(self.handler.mol, path)
            elif ext == "smi":
                from molbuilder.io.smiles_io import write_smiles
                write_smiles(self.handler.mol, path)
            else:
                messagebox.showerror("Error", f"Unsupported format: .{ext}")
                return
            self.status_var.set(f"Saved: {path}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _from_smiles(self):
        dlg = SmilesDialog(self.root)
        if dlg.result:
            smi, name = dlg.result
            try:
                from molbuilder.smiles import parse
                mol = parse(smi)
                mol.name = name or smi
                self.handler.mol = mol
                self.handler.selected_atoms.clear()
                self.status_var.set(f"Built from SMILES: {smi}")
                self._refresh()
            except Exception as e:
                messagebox.showerror("SMILES Error", str(e))

    def _from_formula(self):
        dlg = FormulaDialog(self.root)
        if dlg.result:
            formula, charge = dlg.result
            try:
                from molbuilder.bonding.vsepr import VSEPRMolecule
                vsepr = VSEPRMolecule(formula, charge)
                # Convert VSEPR coords to Molecule
                coords = vsepr.coordinates
                mol = Molecule(formula)
                for sym, pos in coords["atom_positions"]:
                    if sym is not None and pos is not None:
                        mol.add_atom(sym, pos)
                for ai, aj, order in coords["bonds"]:
                    mol.add_bond(ai, aj, order, rotatable=False)
                self.handler.mol = mol
                self.handler.selected_atoms.clear()
                self.status_var.set(f"Built {formula} via VSEPR")
                self._refresh()
            except Exception as e:
                messagebox.showerror("Formula Error", str(e))

    def _preset(self, smiles, name):
        try:
            from molbuilder.smiles import parse
            mol = parse(smiles)
            mol.name = name
            self.handler.mol = mol
            self.handler.selected_atoms.clear()
            self.status_var.set(f"Loaded preset: {name}")
            self._refresh()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _run_analysis(self, name):
        mol = self.handler.mol
        if not mol or len(mol.atoms) == 0:
            self.sidebar.show_results("No molecule loaded.")
            return

        try:
            if name == "Functional Groups":
                from molbuilder.reactions.functional_group_detect import detect_functional_groups
                groups = detect_functional_groups(mol)
                if groups:
                    lines = [f"Found {len(groups)} functional group(s):\n"]
                    for g in groups:
                        lines.append(f"  {g.name} (atom {g.center})")
                    self.sidebar.show_results("\n".join(lines))
                else:
                    self.sidebar.show_results("No functional groups detected.")

            elif name == "Bond Analysis":
                lines = [f"Bond Analysis: {mol.name}\n"]
                lines.append(f"{'Bond':<12} {'Order':>5} {'Length':>8}")
                lines.append("-" * 30)
                for b in mol.bonds:
                    sa = mol.atoms[b.atom_i].symbol
                    sb = mol.atoms[b.atom_j].symbol
                    d = mol.distance(b.atom_i, b.atom_j)
                    sym = {1: "-", 2: "=", 3: "#"}.get(b.order, "?")
                    lines.append(f"{sa}{sym}{sb:<8} {b.order:>5} {d:>7.3f} A")
                self.sidebar.show_results("\n".join(lines))

            elif name == "SMILES":
                from molbuilder.smiles.writer import to_smiles
                smi = to_smiles(mol)
                self.sidebar.show_results(f"SMILES: {smi}")

            elif name == "Stereochemistry":
                lines = ["Stereochemistry Analysis:\n"]
                for i, atom in enumerate(mol.atoms):
                    if atom.symbol == "H":
                        continue
                    nbrs = mol.neighbors(i)
                    if len(nbrs) == 4:
                        desc = mol.assign_rs(i)
                        if desc.name != "NONE":
                            lines.append(f"  {atom.symbol}[{i}]: {desc.name}")
                if len(lines) == 1:
                    lines.append("  No stereocenters found.")
                self.sidebar.show_results("\n".join(lines))

            else:
                self.sidebar.show_results(f"Analysis '{name}' not yet implemented.")

        except Exception as e:
            self.sidebar.show_results(f"Error: {e}")

    # ---- Simulate menu commands ----

    def _sim_md_vibration(self):
        """Launch MD vibration simulation on the current molecule."""
        mol = self.handler.mol
        if not mol or len(mol.atoms) == 0:
            messagebox.showwarning("Simulate", "No molecule loaded.")
            return
        try:
            from molbuilder.visualization.interaction_viz import (
                visualize_md_trajectory, PlaybackConfig,
            )
            self.status_var.set("Running MD simulation...")
            self.root.update()
            config = PlaybackConfig(show_electron_density=False)
            visualize_md_trajectory(mol, n_steps=500, config=config)
            self.status_var.set("MD simulation complete.")
        except Exception as e:
            messagebox.showerror("Simulation Error", str(e))

    def _sim_bond_formation(self):
        """Visualize bond formation between two selected atoms."""
        mol = self.handler.mol
        selected = list(self.handler.selected_atoms)
        if len(selected) != 2:
            messagebox.showinfo(
                "Bond Formation",
                "Select exactly 2 atoms first, then run this command.")
            return
        try:
            from molbuilder.visualization.interaction_viz import (
                visualize_bond_formation, PlaybackConfig,
            )
            self.status_var.set("Simulating bond formation...")
            self.root.update()
            config = PlaybackConfig(show_electron_density=True)
            visualize_bond_formation(
                mol, selected[0], selected[1], config=config)
            self.status_var.set("Bond formation visualization complete.")
        except Exception as e:
            messagebox.showerror("Simulation Error", str(e))

    def _sim_sn2(self):
        """Run an SN2 mechanism demonstration."""
        try:
            from molbuilder.visualization.interaction_viz import (
                visualize_reaction, PlaybackConfig,
            )
            from molbuilder.dynamics.mechanisms import sn2_mechanism
            from molbuilder.molecule.builders import build_ethane

            # Build a simple substrate (methane-like with a "leaving group")
            mol = build_ethane(60.0)
            mechanism = sn2_mechanism(
                substrate_C=0, nucleophile=1, leaving_group=2)
            self.status_var.set("Running SN2 mechanism...")
            self.root.update()
            config = PlaybackConfig(
                show_electron_density=True,
                show_electron_flows=True,
            )
            visualize_reaction(
                mol, mechanism,
                n_steps_per_stage=150, config=config)
            self.status_var.set("SN2 mechanism visualization complete.")
        except Exception as e:
            messagebox.showerror("Simulation Error", str(e))

    def _sim_export(self):
        """Export the current molecule's MD animation to a file."""
        mol = self.handler.mol
        if not mol or len(mol.atoms) == 0:
            messagebox.showwarning("Export", "No molecule loaded.")
            return
        from tkinter import filedialog
        path = filedialog.asksaveasfilename(
            parent=self.root,
            title="Export Animation",
            filetypes=[("GIF", "*.gif"), ("MP4", "*.mp4")],
            defaultextension=".gif",
        )
        if not path:
            return
        try:
            from molbuilder.visualization.interaction_viz import (
                visualize_md_trajectory, PlaybackConfig,
            )
            self.status_var.set(f"Exporting animation to {path}...")
            self.root.update()
            config = PlaybackConfig(
                export_path=path,
                show_electron_density=False,
            )
            viz = visualize_md_trajectory(
                mol, n_steps=300, config=config, show=False)
            viz.export(path)
            self.status_var.set(f"Exported: {path}")
            messagebox.showinfo("Export", f"Animation saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Export Error", str(e))

    def run(self):
        """Start the application main loop."""
        self.root.mainloop()


def launch():
    """Launch the MolBuilder GUI."""
    app = MolBuilderApp()
    app.run()
