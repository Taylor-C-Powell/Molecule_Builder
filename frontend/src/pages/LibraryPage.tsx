import { useState, useEffect, type FormEvent } from "react";
import { useLibrary } from "@/hooks/useLibrary";
import { Button } from "@/components/ui/Button";
import { Input } from "@/components/ui/Input";
import { Card } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import { Alert } from "@/components/ui/Alert";
import type { LibraryMoleculeResponse } from "@/api/types";

function SaveForm({ onSaved }: { onSaved: () => void }) {
  const { save, loading, error } = useLibrary();
  const [smiles, setSmiles] = useState("");
  const [name, setName] = useState("");
  const [tags, setTags] = useState("");
  const [notes, setNotes] = useState("");

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const tagList = tags
      .split(",")
      .map((t) => t.trim())
      .filter(Boolean);
    const result = await save(smiles, name || undefined, tagList.length ? tagList : undefined, notes || undefined);
    if (result) {
      setSmiles("");
      setName("");
      setTags("");
      setNotes("");
      onSaved();
    }
  };

  return (
    <Card>
      <h2 className="text-lg font-semibold mb-4">Save Molecule</h2>
      <form onSubmit={handleSubmit} className="space-y-3">
        <Input label="SMILES" value={smiles} onChange={(e) => setSmiles(e.target.value)} placeholder="e.g. CCO" required />
        <Input label="Name" value={name} onChange={(e) => setName(e.target.value)} placeholder="Optional name" />
        <Input label="Tags" value={tags} onChange={(e) => setTags(e.target.value)} placeholder="Comma-separated tags" />
        <div className="flex flex-col gap-1.5">
          <label className="text-sm font-medium text-text-secondary">Notes</label>
          <textarea
            value={notes}
            onChange={(e) => setNotes(e.target.value)}
            placeholder="Optional notes"
            rows={2}
            className="w-full rounded-[var(--radius-sm)] border border-border bg-bg-card px-3 py-2 text-sm text-text-primary placeholder:text-text-muted focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent transition-colors resize-y"
          />
        </div>
        {error && <Alert variant="error">{error}</Alert>}
        <Button type="submit" disabled={loading || !smiles}>
          {loading ? "Saving..." : "Save to Library"}
        </Button>
      </form>
    </Card>
  );
}

function ImportPanel({ onImported }: { onImported: () => void }) {
  const { importSmiles, loading, error } = useLibrary();
  const [text, setText] = useState("");
  const [tag, setTag] = useState("");
  const [result, setResult] = useState<{ saved: number; duplicates: number; errors: string[] } | null>(null);

  const handleImport = async (e: FormEvent) => {
    e.preventDefault();
    const list = text
      .split("\n")
      .map((s) => s.trim())
      .filter(Boolean);
    if (!list.length) return;
    const res = await importSmiles(list, tag || undefined);
    if (res) {
      setResult(res);
      onImported();
    }
  };

  return (
    <Card>
      <h2 className="text-lg font-semibold mb-4">Bulk Import</h2>
      <form onSubmit={handleImport} className="space-y-3">
        <div className="flex flex-col gap-1.5">
          <label className="text-sm font-medium text-text-secondary">SMILES (one per line)</label>
          <textarea
            value={text}
            onChange={(e) => setText(e.target.value)}
            placeholder={"CCO\nCC(=O)O\nc1ccccc1"}
            rows={4}
            className="w-full rounded-[var(--radius-sm)] border border-border bg-bg-card px-3 py-2 text-sm text-text-primary placeholder:text-text-muted focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent transition-colors resize-y font-mono"
          />
        </div>
        <Input label="Tag (optional)" value={tag} onChange={(e) => setTag(e.target.value)} placeholder="e.g. screening-set" />
        {error && <Alert variant="error">{error}</Alert>}
        {result && (
          <Alert variant="info">
            Saved: {result.saved} | Duplicates: {result.duplicates}
            {result.errors.length > 0 && ` | Errors: ${result.errors.length}`}
          </Alert>
        )}
        <Button type="submit" disabled={loading || !text.trim()}>
          {loading ? "Importing..." : "Import"}
        </Button>
      </form>
    </Card>
  );
}

function MoleculeRow({
  mol,
  onDeleted,
  onUpdated,
}: {
  mol: LibraryMoleculeResponse;
  onDeleted: () => void;
  onUpdated: () => void;
}) {
  const { update, remove, loading } = useLibrary();
  const [expanded, setExpanded] = useState(false);
  const [editing, setEditing] = useState(false);
  const [editName, setEditName] = useState(mol.name ?? "");
  const [editTags, setEditTags] = useState(mol.tags.join(", "));
  const [editNotes, setEditNotes] = useState(mol.notes ?? "");
  const [confirmDelete, setConfirmDelete] = useState(false);

  const handleSave = async () => {
    const tagList = editTags
      .split(",")
      .map((t) => t.trim())
      .filter(Boolean);
    await update(mol.id, { name: editName, tags: tagList, notes: editNotes });
    setEditing(false);
    onUpdated();
  };

  const handleDelete = async () => {
    const ok = await remove(mol.id);
    if (ok) onDeleted();
  };

  return (
    <div className="border border-border rounded-[var(--radius-sm)] bg-bg-card">
      <div
        className="flex items-center justify-between px-4 py-3 cursor-pointer hover:bg-bg-card/80 transition-colors"
        onClick={() => setExpanded(!expanded)}
      >
        <div className="flex items-center gap-3 min-w-0">
          <code className="text-sm font-mono text-accent truncate">{mol.smiles}</code>
          {mol.name && <span className="text-sm text-text-secondary">{mol.name}</span>}
          {mol.tags.map((t) => (
            <Badge key={t} variant="default">{t}</Badge>
          ))}
        </div>
        <span className="text-text-muted text-xs shrink-0">{expanded ? "[-]" : "[+]"}</span>
      </div>

      {expanded && (
        <div className="border-t border-border px-4 py-3 space-y-3">
          {mol.properties && Object.keys(mol.properties).length > 0 && (
            <div className="grid grid-cols-2 sm:grid-cols-4 gap-2 text-xs">
              {Object.entries(mol.properties).map(([k, v]) => (
                <div key={k}>
                  <span className="text-text-muted">{k}: </span>
                  <span className="text-text-primary">{String(v)}</span>
                </div>
              ))}
            </div>
          )}

          {editing ? (
            <div className="space-y-2">
              <Input label="Name" value={editName} onChange={(e) => setEditName(e.target.value)} />
              <Input label="Tags" value={editTags} onChange={(e) => setEditTags(e.target.value)} />
              <div className="flex flex-col gap-1.5">
                <label className="text-sm font-medium text-text-secondary">Notes</label>
                <textarea
                  value={editNotes}
                  onChange={(e) => setEditNotes(e.target.value)}
                  rows={2}
                  className="w-full rounded-[var(--radius-sm)] border border-border bg-bg-card px-3 py-2 text-sm text-text-primary focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent transition-colors resize-y"
                />
              </div>
              <div className="flex gap-2">
                <Button size="sm" onClick={handleSave} disabled={loading}>Save</Button>
                <Button size="sm" variant="ghost" onClick={() => setEditing(false)}>Cancel</Button>
              </div>
            </div>
          ) : (
            <div className="flex gap-2">
              <Button size="sm" variant="ghost" onClick={() => setEditing(true)}>Edit</Button>
              {confirmDelete ? (
                <>
                  <Button size="sm" variant="ghost" onClick={handleDelete} disabled={loading} className="text-red">
                    Confirm Delete
                  </Button>
                  <Button size="sm" variant="ghost" onClick={() => setConfirmDelete(false)}>Cancel</Button>
                </>
              ) : (
                <Button size="sm" variant="ghost" onClick={() => setConfirmDelete(true)}>Delete</Button>
              )}
            </div>
          )}

          {mol.notes && !editing && <p className="text-xs text-text-muted">{mol.notes}</p>}
        </div>
      )}
    </div>
  );
}

export default function LibraryPage() {
  const { fetchList, loading, error, list } = useLibrary();
  const [tagFilter, setTagFilter] = useState("");
  const [search, setSearch] = useState("");
  const [page, setPage] = useState(1);

  const reload = () => {
    fetchList({
      tag: tagFilter || undefined,
      search: search || undefined,
      page,
    });
  };

  useEffect(() => {
    reload();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [page]);

  const handleSearch = (e: FormEvent) => {
    e.preventDefault();
    setPage(1);
    reload();
  };

  const totalPages = list ? Math.ceil(list.total / list.per_page) : 0;

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-2xl font-bold mb-1">Molecule Library</h1>
        <p className="text-sm text-text-secondary">
          Save, organize, and manage your molecules.
        </p>
      </div>

      <div className="grid md:grid-cols-2 gap-6">
        <SaveForm onSaved={reload} />
        <ImportPanel onImported={reload} />
      </div>

      <Card>
        <form onSubmit={handleSearch} className="flex flex-wrap gap-3 mb-4">
          <Input
            placeholder="Search SMILES or name..."
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            className="flex-1 min-w-[200px]"
          />
          <Input
            placeholder="Filter by tag"
            value={tagFilter}
            onChange={(e) => setTagFilter(e.target.value)}
            className="w-40"
          />
          <Button type="submit" disabled={loading}>
            {loading ? "Loading..." : "Search"}
          </Button>
        </form>

        {error && <Alert variant="error" className="mb-4">{error}</Alert>}

        {list && (
          <>
            <p className="text-xs text-text-muted mb-3">
              {list.total} molecule{list.total !== 1 ? "s" : ""} found
            </p>
            <div className="space-y-2">
              {list.molecules.map((mol) => (
                <MoleculeRow key={mol.id} mol={mol} onDeleted={reload} onUpdated={reload} />
              ))}
            </div>

            {totalPages > 1 && (
              <div className="flex items-center justify-center gap-2 mt-4">
                <Button size="sm" variant="ghost" disabled={page <= 1} onClick={() => setPage(page - 1)}>
                  Prev
                </Button>
                <span className="text-sm text-text-secondary">
                  Page {page} of {totalPages}
                </span>
                <Button size="sm" variant="ghost" disabled={page >= totalPages} onClick={() => setPage(page + 1)}>
                  Next
                </Button>
              </div>
            )}
          </>
        )}
      </Card>
    </div>
  );
}
