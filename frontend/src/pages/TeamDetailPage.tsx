import { useEffect, useState } from "react";
import { useParams, Link, useNavigate } from "react-router-dom";
import { useTeams } from "@/hooks/useTeams";
import { MembersPanel } from "@/components/teams/MembersPanel";
import { TeamLibraryPanel } from "@/components/teams/TeamLibraryPanel";
import { Button } from "@/components/ui/Button";
import { Input } from "@/components/ui/Input";
import { Alert } from "@/components/ui/Alert";
import { cn } from "@/lib/cn";

type Tab = "members" | "library";

export default function TeamDetailPage() {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const teamId = Number(id);
  const { fetchTeam, updateTeam, deleteTeam, selectedTeam, loading, error } =
    useTeams();

  const [tab, setTab] = useState<Tab>("members");
  const [editing, setEditing] = useState(false);
  const [editName, setEditName] = useState("");
  const [confirmDelete, setConfirmDelete] = useState(false);

  useEffect(() => {
    if (teamId) fetchTeam(teamId);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [teamId]);

  useEffect(() => {
    if (selectedTeam) setEditName(selectedTeam.name);
  }, [selectedTeam]);

  const handleUpdate = async () => {
    if (!editName) return;
    const result = await updateTeam(teamId, editName);
    if (result) setEditing(false);
  };

  const handleDelete = async () => {
    const ok = await deleteTeam(teamId);
    if (ok) navigate("/teams");
  };

  if (!selectedTeam && !loading) {
    return (
      <div className="space-y-4">
        <Link
          to="/teams"
          className="text-sm text-accent hover:underline"
        >
          Back to Teams
        </Link>
        {error && <Alert variant="error">{error}</Alert>}
        <p className="text-text-muted">Team not found.</p>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <Link
        to="/teams"
        className="text-sm text-accent hover:underline"
      >
        Back to Teams
      </Link>

      {selectedTeam && (
        <div className="flex items-start justify-between">
          <div>
            {editing ? (
              <div className="flex items-center gap-2">
                <Input
                  value={editName}
                  onChange={(e) => setEditName(e.target.value)}
                  className="text-lg font-bold"
                />
                <Button size="sm" onClick={handleUpdate} disabled={loading}>
                  Save
                </Button>
                <Button
                  size="sm"
                  variant="ghost"
                  onClick={() => {
                    setEditing(false);
                    setEditName(selectedTeam.name);
                  }}
                >
                  Cancel
                </Button>
              </div>
            ) : (
              <div className="flex items-center gap-3">
                <h1 className="text-2xl font-bold">{selectedTeam.name}</h1>
                <Button size="sm" variant="ghost" onClick={() => setEditing(true)}>
                  Edit
                </Button>
              </div>
            )}
            <p className="text-xs text-text-muted font-mono mt-1">
              {selectedTeam.slug}
            </p>
            <p className="text-xs text-text-secondary mt-1">
              Owner: {selectedTeam.owner_email}
            </p>
          </div>

          <div>
            {confirmDelete ? (
              <div className="flex gap-2">
                <Button
                  size="sm"
                  variant="ghost"
                  onClick={handleDelete}
                  disabled={loading}
                  className="text-red"
                >
                  Confirm Delete
                </Button>
                <Button
                  size="sm"
                  variant="ghost"
                  onClick={() => setConfirmDelete(false)}
                >
                  Cancel
                </Button>
              </div>
            ) : (
              <Button
                size="sm"
                variant="ghost"
                onClick={() => setConfirmDelete(true)}
              >
                Delete Team
              </Button>
            )}
          </div>
        </div>
      )}

      {error && <Alert variant="error">{error}</Alert>}

      {/* Tabs */}
      <div className="flex border-b border-border">
        {(["members", "library"] as const).map((t) => (
          <button
            key={t}
            onClick={() => setTab(t)}
            className={cn(
              "px-4 py-2 text-sm font-medium border-b-2 transition-colors bg-transparent cursor-pointer",
              tab === t
                ? "border-accent text-text-primary"
                : "border-transparent text-text-secondary hover:text-text-primary",
            )}
          >
            {t === "members" ? "Members" : "Library"}
          </button>
        ))}
      </div>

      {selectedTeam && tab === "members" && (
        <MembersPanel
          teamId={teamId}
          ownerEmail={selectedTeam.owner_email}
        />
      )}

      {selectedTeam && tab === "library" && (
        <TeamLibraryPanel teamId={teamId} />
      )}
    </div>
  );
}
