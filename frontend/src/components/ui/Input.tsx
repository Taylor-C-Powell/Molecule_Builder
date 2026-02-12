import { cn } from "@/lib/cn";
import type { InputHTMLAttributes } from "react";

interface InputProps extends InputHTMLAttributes<HTMLInputElement> {
  label?: string;
  error?: string;
}

export function Input({ label, error, className, id, ...props }: InputProps) {
  const inputId = id ?? label?.toLowerCase().replace(/\s+/g, "-");
  return (
    <div className="flex flex-col gap-1.5">
      {label && (
        <label htmlFor={inputId} className="text-sm font-medium text-text-secondary">
          {label}
        </label>
      )}
      <input
        id={inputId}
        className={cn(
          "w-full rounded-[var(--radius-sm)] border border-border bg-bg-card px-3 py-2 text-sm text-text-primary",
          "placeholder:text-text-muted",
          "focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent",
          "transition-colors",
          error && "border-red",
          className,
        )}
        {...props}
      />
      {error && <p className="text-xs text-red">{error}</p>}
    </div>
  );
}
