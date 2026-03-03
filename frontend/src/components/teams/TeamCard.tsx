import { Link } from "react-router-dom";
import { Card } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import type { TeamResponse } from "@/api/types";

export function TeamCard({ team, currentEmail }: { team: TeamResponse; currentEmail?: string }) {
  const isOwner = currentEmail === team.owner_email;

  return (
    <Link to={`/teams/${team.id}`} className="block no-underline">
      <Card className="hover:border-accent transition-colors cursor-pointer">
        <div className="flex items-center justify-between">
          <div>
            <h3 className="text-base font-semibold text-text-primary">{team.name}</h3>
            <p className="text-xs text-text-muted font-mono">{team.slug}</p>
          </div>
          <div className="flex items-center gap-2">
            {isOwner && <Badge variant="default">Owner</Badge>}
          </div>
        </div>
        <p className="text-xs text-text-secondary mt-2">
          Owner: {team.owner_email}
        </p>
      </Card>
    </Link>
  );
}
