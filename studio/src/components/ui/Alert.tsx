import { cn } from "@/lib/cn";
import type { HTMLAttributes } from "react";

type Variant = "info" | "success" | "warning" | "error";

interface AlertProps extends HTMLAttributes<HTMLDivElement> {
  variant?: Variant;
}

const variantStyles: Record<Variant, string> = {
  info: "border-accent/30 bg-accent/5 text-text-secondary",
  success: "border-green/30 bg-green/5 text-green",
  warning: "border-yellow/30 bg-yellow/5 text-yellow",
  error: "border-red/30 bg-red/5 text-red",
};

export function Alert({ variant = "info", className, ...props }: AlertProps) {
  return (
    <div
      className={cn(
        "px-4 py-3 text-sm border rounded-[var(--radius-sm)]",
        variantStyles[variant],
        className,
      )}
      {...props}
    />
  );
}
