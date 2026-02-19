import { cn } from "@/lib/cn";
import type { InputHTMLAttributes } from "react";

interface InputProps extends InputHTMLAttributes<HTMLInputElement> {
  label?: string;
  error?: string;
}

export function Input({ label, error, className, ...props }: InputProps) {
  return (
    <div className="space-y-1">
      {label && (
        <label className="block text-xs font-medium text-text-secondary">
          {label}
        </label>
      )}
      <input
        className={cn(
          "w-full px-3 py-2 text-sm bg-bg-card text-text-primary border border-border",
          "rounded-[var(--radius-sm)] outline-none transition-colors",
          "focus:border-accent placeholder:text-text-muted",
          error && "border-red",
          className,
        )}
        {...props}
      />
      {error && <p className="text-xs text-red">{error}</p>}
    </div>
  );
}
