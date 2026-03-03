import { useEffect } from "react";
import { useTeams } from "@/hooks/useTeams";
import { useAuthStore } from "@/stores/auth-store";
import { CreateTeamForm } from "@/components/teams/CreateTeamForm";
import { TeamCard } from "@/components/teams/TeamCard";
import { Alert } from "@/components/ui/Alert";

export default function TeamsPage() {
  const { fetchTeams, createTeam, teams, loading, error } = useTeams();
  const email = useAuthStore((s) => s.email);

  useEffect(() => {
    fetchTeams();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-2xl font-bold mb-1">Teams</h1>
        <p className="text-sm text-text-secondary">
          Create and manage teams for collaborative molecule work.
        </p>
      </div>

      <CreateTeamForm
        onCreated={() => fetchTeams()}
        loading={loading}
        error={error}
        onCreate={createTeam}
      />

      {error && <Alert variant="error">{error}</Alert>}

      {teams.length > 0 ? (
        <div className="space-y-3">
          {teams.map((team) => (
            <TeamCard key={team.id} team={team} currentEmail={email ?? undefined} />
          ))}
        </div>
      ) : (
        !loading && (
          <p className="text-sm text-text-muted">
            No teams yet. Create one to get started.
          </p>
        )
      )}
    </div>
  );
}
