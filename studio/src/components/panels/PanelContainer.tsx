import { useState, type ReactNode } from "react";
import { cn } from "@/lib/cn";

interface PanelContainerProps {
  title: string;
  defaultOpen?: boolean;
  children: ReactNode;
}

export function PanelContainer({
  title,
  defaultOpen = true,
  children,
}: PanelContainerProps) {
  const [open, setOpen] = useState(defaultOpen);

  return (
    <div className="border-b border-border">
      <button
        className="w-full flex items-center justify-between px-3 py-2 text-xs font-semibold text-text-secondary hover:text-text-primary transition-colors cursor-pointer"
        onClick={() => setOpen(!open)}
      >
        {title}
        <svg
          className={cn("w-3 h-3 transition-transform", open && "rotate-180")}
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
          strokeWidth={2}
        >
          <path strokeLinecap="round" strokeLinejoin="round" d="M19 9l-7 7-7-7" />
        </svg>
      </button>
      {open && <div className="px-3 pb-3">{children}</div>}
    </div>
  );
}
