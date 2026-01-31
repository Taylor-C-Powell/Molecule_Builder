"""Interactive menu system for Molecule Builder.

Provides a CLI with numbered options for each demo module,
an 'all' mode, and command-line argument support.

Migrated from legacy/main.py.
"""

import sys

from molbuilder.cli.demos import (
    demo_bohr_model,
    demo_quantum_model,
    demo_element_data,
    demo_lewis_structure,
    demo_vsepr_model,
    demo_covalent_bonds,
    demo_molecular_conformations,
    demo_amino_acids,
    demo_visualization,
    demo_quantum_visualization,
)
from molbuilder.cli.wizard import wizard_main


# ===================================================================
# Menu
# ===================================================================

DEMOS = [
    ("Bohr Atomic Model",                  demo_bohr_model),
    ("Quantum Mechanical Atom",            demo_quantum_model),
    ("Element Data",                        demo_element_data),
    ("Lewis Structures",                    demo_lewis_structure),
    ("VSEPR Molecular Geometry",            demo_vsepr_model),
    ("Covalent Bonds",                      demo_covalent_bonds),
    ("Molecular Conformations",             demo_molecular_conformations),
    ("Amino Acids & Functional Groups",     demo_amino_acids),
    ("3D Molecule Visualization",           demo_visualization),
    ("Quantum Orbital Visualization",       demo_quantum_visualization),
    ("Molecule Builder Wizard",             wizard_main),
]


def print_banner():
    print()
    print("=" * 60)
    print("         MOLECULE BUILDER  --  Alpha 1.0")
    print("         Taylor C. Powell")
    print("=" * 60)
    print()


def print_menu():
    print("  Available modules:")
    print()
    for i, (name, _) in enumerate(DEMOS, 1):
        print(f"    [{i:>2}] {name}")
    print()
    print(f"    [ a] Run all text demos (1-8)")
    print(f"    [ q] Quit")
    print()


def run_all_text():
    """Run all non-visualization demos in sequence."""
    for name, func in DEMOS[:8]:
        try:
            func()
        except Exception as e:
            print(f"  ERROR in {name}: {e}")
            print()


def main():
    print_banner()

    if len(sys.argv) > 1:
        arg = sys.argv[1].lower()
        if arg == "all":
            run_all_text()
            return
        if arg == "q":
            return
        try:
            idx = int(arg) - 1
            if 0 <= idx < len(DEMOS):
                DEMOS[idx][1]()
                return
        except ValueError:
            pass
        print(f"  Unknown argument: {sys.argv[1]}")
        print(f"  Usage: python main.py [1-{len(DEMOS)} | all | q]")
        return

    # Interactive menu loop
    while True:
        print_menu()
        choice = input("  Select option: ").strip().lower()

        if choice == "q":
            print("  Goodbye!")
            break

        if choice == "a":
            run_all_text()
            continue

        try:
            idx = int(choice) - 1
            if 0 <= idx < len(DEMOS):
                print()
                try:
                    DEMOS[idx][1]()
                except Exception as e:
                    print(f"  ERROR: {e}")
                    import traceback
                    traceback.print_exc()
                    print()
            else:
                print(f"  Please enter 1-{len(DEMOS)}, 'a', or 'q'.")
        except ValueError:
            print(f"  Please enter 1-{len(DEMOS)}, 'a', or 'q'.")


if __name__ == "__main__":
    main()
