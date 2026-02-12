import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";

interface FunctionalGroupsProps {
  groups: string[];
}

export function FunctionalGroups({ groups }: FunctionalGroupsProps) {
  if (!groups.length) return null;

  return (
    <Card>
      <CardHeader>
        <CardTitle>Functional Groups</CardTitle>
      </CardHeader>
      <div className="flex flex-wrap gap-2">
        {groups.map((g) => (
          <Badge key={g} variant="accent">
            {g}
          </Badge>
        ))}
      </div>
    </Card>
  );
}
