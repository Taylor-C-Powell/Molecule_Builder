import { useState, useEffect, type FormEvent } from "react";
import { Button } from "@/components/ui/Button";
import { Input } from "@/components/ui/Input";
import { Card } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import { Alert } from "@/components/ui/Alert";
import { useTeams } from "@/hooks/useTeams";
import type { MemberResponse } from "@/api/types";

const ROLES = ["member", "admin"] as const;

export function MembersPanel({
  teamId,
  ownerEmail: _ownerEmail,
}: {
  teamId: number;
  ownerEmail: string;
}) {
  const {
    fetchMembers,
    addMember,
    updateMemberRole,
    removeMember,
    members,
    loading,
    error,
  } = useTeams();

  const [newEmail, setNewEmail] = useState("");
  const [newRole, setNewRole] = useState<string>("member");
  const [confirmRemove, setConfirmRemove] = useState<string | null>(null);

  useEffect(() => {
    fetchMembers(teamId);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [teamId]);

  const handleAdd = async (e: FormEvent) => {
    e.preventDefault();
    if (!newEmail) return;
    const result = await addMember(teamId, newEmail, newRole);
    if (result) {
      setNewEmail("");
      setNewRole("member");
      fetchMembers(teamId);
    }
  };

  const handleRoleChange = async (member: MemberResponse, role: string) => {
    await updateMemberRole(teamId, member.user_email, role);
    fetchMembers(teamId);
  };

  const handleRemove = async (email: string) => {
    const ok = await removeMember(teamId, email);
    if (ok) {
      setConfirmRemove(null);
      fetchMembers(teamId);
    }
  };

  return (
    <Card>
      <h3 className="text-lg font-semibold mb-4">Members</h3>

      <form onSubmit={handleAdd} className="flex flex-wrap gap-3 mb-4">
        <Input
          placeholder="Email address"
          value={newEmail}
          onChange={(e) => setNewEmail(e.target.value)}
          className="flex-1 min-w-[200px]"
          required
        />
        <select
          value={newRole}
          onChange={(e) => setNewRole(e.target.value)}
          className="rounded-[var(--radius-sm)] border border-border bg-bg-card px-3 py-2 text-sm text-text-primary focus:outline-none focus:border-accent"
        >
          {ROLES.map((r) => (
            <option key={r} value={r}>
              {r}
            </option>
          ))}
        </select>
        <Button type="submit" disabled={loading || !newEmail}>
          {loading ? "Adding..." : "Add Member"}
        </Button>
      </form>

      {error && <Alert variant="error" className="mb-4">{error}</Alert>}

      <div className="space-y-2">
        {members.map((m) => {
          const isOwner = m.team_role === "owner";
          return (
            <div
              key={m.id}
              className="flex items-center justify-between border border-border rounded-[var(--radius-sm)] px-4 py-3"
            >
              <div className="flex items-center gap-3">
                <span className="text-sm text-text-primary">{m.user_email}</span>
                <Badge variant={isOwner ? "default" : "default"}>
                  {m.team_role}
                </Badge>
              </div>
              <div className="flex items-center gap-2">
                {!isOwner && (
                  <>
                    <select
                      value={m.team_role}
                      onChange={(e) => handleRoleChange(m, e.target.value)}
                      disabled={loading}
                      className="rounded-[var(--radius-sm)] border border-border bg-bg-card px-2 py-1 text-xs text-text-primary focus:outline-none focus:border-accent"
                    >
                      {ROLES.map((r) => (
                        <option key={r} value={r}>
                          {r}
                        </option>
                      ))}
                    </select>
                    {confirmRemove === m.user_email ? (
                      <div className="flex gap-1">
                        <Button
                          size="sm"
                          variant="ghost"
                          onClick={() => handleRemove(m.user_email)}
                          disabled={loading}
                          className="text-red"
                        >
                          Confirm
                        </Button>
                        <Button
                          size="sm"
                          variant="ghost"
                          onClick={() => setConfirmRemove(null)}
                        >
                          Cancel
                        </Button>
                      </div>
                    ) : (
                      <Button
                        size="sm"
                        variant="ghost"
                        onClick={() => setConfirmRemove(m.user_email)}
                      >
                        Remove
                      </Button>
                    )}
                  </>
                )}
              </div>
            </div>
          );
        })}
        {members.length === 0 && !loading && (
          <p className="text-sm text-text-muted">No members yet.</p>
        )}
      </div>
    </Card>
  );
}
