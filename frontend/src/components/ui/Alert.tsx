import { cn } from "@/lib/cn";

type AlertVariant = "info" | "success" | "warning" | "error";

const variantStyles: Record<AlertVariant, string> = {
  info: "border-accent/30 bg-accent/5 text-accent",
  success: "border-green/30 bg-green/5 text-green",
  warning: "border-yellow/30 bg-yellow/5 text-yellow",
  error: "border-red/30 bg-red/5 text-red",
};

interface AlertProps {
  variant?: AlertVariant;
  className?: string;
  children: React.ReactNode;
}

export function Alert({ variant = "info", className, children }: AlertProps) {
  return (
    <div
      role="alert"
      className={cn(
        "rounded-[var(--radius-sm)] border px-4 py-3 text-sm",
        variantStyles[variant],
        className,
      )}
    >
      {children}
    </div>
  );
}
