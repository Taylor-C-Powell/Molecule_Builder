export const API_BASE =
  import.meta.env.VITE_API_URL ??
  (import.meta.env.DEV
    ? "/api"
    : "https://molbuilder-api-production.up.railway.app/api");

export const APP_NAME = "MolBuilder";

export const TIERS = {
  free: { label: "Free", rpm: 10, expensivePerHour: 5, price: 0 },
  pro: { label: "Pro", rpm: 60, expensivePerHour: 30, price: 49 },
  team: { label: "Team", rpm: 120, expensivePerHour: 60, price: 199 },
  academic: { label: "Academic", rpm: 30, expensivePerHour: 15, price: 0 },
  enterprise: { label: "Enterprise", rpm: 300, expensivePerHour: 200, price: 499 },
} as const;

export type TierName = keyof typeof TIERS;
