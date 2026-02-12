import { cn } from "@/lib/cn";

type BadgeVariant = "default" | "success" | "warning" | "danger" | "accent";

const variants: Record<BadgeVariant, string> = {
  default: "bg-border text-text-secondary",
  success: "bg-green/15 text-green",
  warning: "bg-yellow/15 text-yellow",
  danger: "bg-red/15 text-red",
  accent: "bg-accent/15 text-accent",
};

interface BadgeProps {
  variant?: BadgeVariant;
  className?: string;
  children: React.ReactNode;
}

export function Badge({ variant = "default", className, children }: BadgeProps) {
  return (
    <span
      className={cn(
        "inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium",
        variants[variant],
        className,
      )}
    >
      {children}
    </span>
  );
}
