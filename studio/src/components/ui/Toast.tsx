import { useToastStore, type ToastType } from "@/stores/toast-store";
import { cn } from "@/lib/cn";

const TYPE_STYLES: Record<ToastType, string> = {
  info: "bg-blue-900/80 border-blue-500/40 text-blue-200",
  success: "bg-green-900/80 border-green-500/40 text-green-200",
  error: "bg-red-900/80 border-red-500/40 text-red-200",
  warning: "bg-amber-900/80 border-amber-500/40 text-amber-200",
};

export function ToastContainer() {
  const toasts = useToastStore((s) => s.toasts);
  const removeToast = useToastStore((s) => s.removeToast);

  if (toasts.length === 0) return null;

  return (
    <div className="fixed bottom-10 left-1/2 -translate-x-1/2 z-[100] flex flex-col gap-1.5 pointer-events-none">
      {toasts.map((toast) => (
        <div
          key={toast.id}
          className={cn(
            "px-4 py-2 rounded-[var(--radius-sm)] border text-xs font-mono shadow-lg pointer-events-auto cursor-pointer animate-[fadeIn_0.2s_ease-out]",
            TYPE_STYLES[toast.type],
          )}
          onClick={() => removeToast(toast.id)}
        >
          {toast.message}
        </div>
      ))}
    </div>
  );
}
