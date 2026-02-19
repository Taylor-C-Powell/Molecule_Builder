import { cn } from "@/lib/cn";
import type { HTMLAttributes } from "react";

type Variant = "default" | "success" | "warning" | "danger" | "accent";

interface BadgeProps extends HTMLAttributes<HTMLSpanElement> {
  variant?: Variant;
}

const variantStyles: Record<Variant, string> = {
  default: "bg-border text-text-secondary",
  success: "bg-green/15 text-green",
  warning: "bg-yellow/15 text-yellow",
  danger: "bg-red/15 text-red",
  accent: "bg-accent/15 text-accent",
};

export function Badge({ variant = "default", className, ...props }: BadgeProps) {
  return (
    <span
      className={cn(
        "inline-flex items-center px-2 py-0.5 text-xs font-medium rounded-full",
        variantStyles[variant],
        className,
      )}
      {...props}
    />
  );
}
