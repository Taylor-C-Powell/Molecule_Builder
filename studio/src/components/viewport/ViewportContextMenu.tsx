import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import { ContextMenu, type ContextMenuItem } from "./ContextMenu";
import { PLACEABLE_ELEMENTS } from "@/types/editor";

export function ViewportContextMenu() {
  const contextMenu = useEditorStore((s) => s.contextMenu);
  const closeContextMenu = useEditorStore((s) => s.closeContextMenu);

  if (!contextMenu) return null;

  const entry = useWorkspaceStore.getState().getActive();
  if (!entry?.structure) return null;

  const items: ContextMenuItem[] = [];

  if (contextMenu.type === "atom") {
    const atom = entry.structure.atoms[contextMenu.index];
    if (!atom) return null;

    items.push({
      label: `${atom.symbol}#${atom.index}`,
      action: () => {},
      disabled: true,
    });
    items.push({ label: "", action: () => {}, separator: true });

    // Select
    items.push({
      label: "Select",
      action: () => useEditorStore.getState().selectAtom(contextMenu.index),
    });

    // Change element submenu -- pick top 4 different elements
    const alternates = PLACEABLE_ELEMENTS.filter((el) => el !== atom.symbol).slice(0, 4);
    for (const el of alternates) {
      items.push({
        label: `Change to ${el}`,
        action: () => {
          useHistoryStore.getState().pushSnapshot(entry.structure!, `Change to ${el}`);
          useWorkspaceStore.getState().replaceAtomElement(contextMenu.index, el);
        },
      });
    }

    items.push({ label: "", action: () => {}, separator: true });

    // Select connected component
    items.push({
      label: "Select connected",
      action: () => {
        const structure = entry.structure!;
        const visited = new Set<number>();
        const queue = [contextMenu.index];
        while (queue.length > 0) {
          const idx = queue.pop()!;
          if (visited.has(idx)) continue;
          visited.add(idx);
          for (const b of structure.bonds) {
            if (b.atom_i === idx && !visited.has(b.atom_j)) queue.push(b.atom_j);
            if (b.atom_j === idx && !visited.has(b.atom_i)) queue.push(b.atom_i);
          }
        }
        const editor = useEditorStore.getState();
        for (const idx of visited) {
          editor.selectAtom(idx, true);
        }
      },
    });

    // Copy position
    items.push({
      label: "Copy position",
      action: () => {
        const pos = atom.position;
        navigator.clipboard.writeText(`${pos[0].toFixed(4)}, ${pos[1].toFixed(4)}, ${pos[2].toFixed(4)}`);
      },
    });

    items.push({ label: "", action: () => {}, separator: true });

    items.push({
      label: "Delete atom",
      action: () => {
        useHistoryStore.getState().pushSnapshot(entry.structure!, "Delete atom");
        useWorkspaceStore.getState().removeAtoms(new Set([contextMenu.index]));
        useEditorStore.getState().clearSelection();
      },
    });
  } else if (contextMenu.type === "bond" && contextMenu.bondKey) {
    const bond = entry.structure.bonds.find(
      (b) => `${b.atom_i}-${b.atom_j}` === contextMenu.bondKey || `${b.atom_j}-${b.atom_i}` === contextMenu.bondKey,
    );
    if (!bond) return null;

    const atomA = entry.structure.atoms[bond.atom_i];
    const atomB = entry.structure.atoms[bond.atom_j];
    const orderLabels = ["", "Single", "Double", "Triple"];

    items.push({
      label: `${atomA?.symbol}#${bond.atom_i} - ${atomB?.symbol}#${bond.atom_j} (${orderLabels[bond.order]})`,
      action: () => {},
      disabled: true,
    });
    items.push({ label: "", action: () => {}, separator: true });

    items.push({
      label: "Select bond",
      action: () => useEditorStore.getState().selectBond(contextMenu.bondKey!),
    });

    // Bond order changes
    for (const order of [1, 2, 3] as const) {
      if (order !== bond.order) {
        items.push({
          label: `Set ${(orderLabels[order] ?? "").toLowerCase()} bond`,
          action: () => {
            useHistoryStore.getState().pushSnapshot(entry.structure!, `Set bond order ${order}`);
            // Cycle to target -- just set directly via workspace
            const structure = entry.structure!;
            const bondIdx = structure.bonds.indexOf(bond);
            if (bondIdx >= 0) {
              const newBonds = structure.bonds.map((b, i) =>
                i === bondIdx ? { ...b, order, rotatable: order === 1 } : b,
              );
              useWorkspaceStore.getState().updateActiveStructure({ ...structure, bonds: newBonds });
            }
          },
        });
      }
    }

    items.push({ label: "", action: () => {}, separator: true });

    items.push({
      label: "Select both atoms",
      action: () => {
        const editor = useEditorStore.getState();
        editor.selectAtom(bond.atom_i);
        editor.selectAtom(bond.atom_j, true);
      },
    });

    items.push({
      label: "Delete bond",
      action: () => {
        useHistoryStore.getState().pushSnapshot(entry.structure!, "Delete bond");
        useWorkspaceStore.getState().removeBonds(new Set([contextMenu.bondKey!]));
        useEditorStore.getState().clearSelection();
      },
    });
  }

  return <ContextMenu x={contextMenu.x} y={contextMenu.y} items={items} onClose={closeContextMenu} />;
}
