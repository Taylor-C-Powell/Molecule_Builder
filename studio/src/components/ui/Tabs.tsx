import { createContext, useContext, useState, useCallback } from "react";
import { cn } from "@/lib/cn";
import type { HTMLAttributes, ReactNode } from "react";

interface TabsContextValue {
  active: string;
  setActive: (id: string) => void;
}

const TabsContext = createContext<TabsContextValue | null>(null);

function useTabs() {
  const ctx = useContext(TabsContext);
  if (!ctx) throw new Error("Tabs components must be used within <Tabs>");
  return ctx;
}

interface TabsProps extends HTMLAttributes<HTMLDivElement> {
  defaultValue: string;
}

export function Tabs({ defaultValue, children, ...props }: TabsProps) {
  const [active, setActive] = useState(defaultValue);
  return (
    <TabsContext.Provider value={{ active, setActive }}>
      <div {...props}>{children}</div>
    </TabsContext.Provider>
  );
}

export function TabsList({ className, ...props }: HTMLAttributes<HTMLDivElement>) {
  return (
    <div
      className={cn(
        "flex gap-1 border-b border-border px-2",
        className,
      )}
      {...props}
    />
  );
}

interface TabsTriggerProps extends HTMLAttributes<HTMLButtonElement> {
  value: string;
  children: ReactNode;
}

export function TabsTrigger({ value, className, ...props }: TabsTriggerProps) {
  const { active, setActive } = useTabs();
  const isActive = active === value;

  const handleClick = useCallback(() => setActive(value), [setActive, value]);

  return (
    <button
      className={cn(
        "px-3 py-2 text-xs font-medium transition-colors cursor-pointer",
        isActive
          ? "text-accent border-b-2 border-accent -mb-px"
          : "text-text-muted hover:text-text-secondary",
        className,
      )}
      onClick={handleClick}
      {...props}
    />
  );
}

interface TabsContentProps extends HTMLAttributes<HTMLDivElement> {
  value: string;
}

export function TabsContent({ value, ...props }: TabsContentProps) {
  const { active } = useTabs();
  if (active !== value) return null;
  return <div {...props} />;
}
