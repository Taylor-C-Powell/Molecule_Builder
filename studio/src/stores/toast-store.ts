import { create } from "zustand";

export type ToastType = "info" | "success" | "error" | "warning";

export interface Toast {
  id: number;
  message: string;
  type: ToastType;
}

interface ToastState {
  toasts: Toast[];
  addToast: (message: string, type?: ToastType) => void;
  removeToast: (id: number) => void;
}

let nextId = 0;

export const useToastStore = create<ToastState>((set, get) => ({
  toasts: [],

  addToast(message, type = "info") {
    const id = nextId++;
    const toast: Toast = { id, message, type };
    set({ toasts: [...get().toasts, toast] });

    // Auto-dismiss after 3 seconds
    setTimeout(() => {
      const current = get().toasts;
      set({ toasts: current.filter((t) => t.id !== id) });
    }, 3000);
  },

  removeToast(id) {
    set({ toasts: get().toasts.filter((t) => t.id !== id) });
  },
}));
