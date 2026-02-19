export const API_BASE =
  import.meta.env.VITE_API_URL ??
  (import.meta.env.DEV
    ? "/api"
    : "https://molbuilder-api-production.up.railway.app/api");

export const APP_NAME = "MolBuilder Studio";
