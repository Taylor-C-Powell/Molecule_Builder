import { create } from "zustand";
import { getToken } from "@/api/auth";

interface AuthState {
  apiKey: string | null;
  jwt: string | null;
  jwtExpiresAt: number | null;
  email: string | null;
  tier: string | null;
  isAuthenticated: boolean;

  setApiKey: (key: string, email: string, tier: string) => void;
  login: (apiKey: string) => Promise<void>;
  logout: () => void;
  refreshToken: () => Promise<void>;
  getValidToken: () => Promise<string | null>;
}

const STORAGE_KEY = "molbuilder_auth";

function loadPersistedState() {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (raw) {
      const parsed = JSON.parse(raw);
      return {
        apiKey: parsed.apiKey ?? null,
        email: parsed.email ?? null,
        tier: parsed.tier ?? null,
      };
    }
  } catch {
    // corrupt storage
  }
  return { apiKey: null, email: null, tier: null };
}

function persist(apiKey: string | null, email: string | null, tier: string | null) {
  if (apiKey) {
    localStorage.setItem(STORAGE_KEY, JSON.stringify({ apiKey, email, tier }));
  } else {
    localStorage.removeItem(STORAGE_KEY);
  }
}

const initial = loadPersistedState();

export const useAuthStore = create<AuthState>((set, get) => ({
  apiKey: initial.apiKey,
  jwt: null,
  jwtExpiresAt: null,
  email: initial.email,
  tier: initial.tier,
  isAuthenticated: !!initial.apiKey,

  setApiKey(key, email, tier) {
    persist(key, email, tier);
    set({ apiKey: key, email, tier, isAuthenticated: true });
  },

  async login(apiKey) {
    const res = await getToken(apiKey);
    const expiresAt = Date.now() + res.expires_in * 1000;
    set({ jwt: res.access_token, jwtExpiresAt: expiresAt });
  },

  logout() {
    persist(null, null, null);
    set({
      apiKey: null,
      jwt: null,
      jwtExpiresAt: null,
      email: null,
      tier: null,
      isAuthenticated: false,
    });
  },

  async refreshToken() {
    const { apiKey } = get();
    if (!apiKey) return;
    try {
      const res = await getToken(apiKey);
      const expiresAt = Date.now() + res.expires_in * 1000;
      set({ jwt: res.access_token, jwtExpiresAt: expiresAt });
    } catch {
      // API key might be revoked
      get().logout();
    }
  },

  async getValidToken() {
    const { jwt, jwtExpiresAt, apiKey } = get();
    // Refresh if expiring within 5 minutes
    if (jwt && jwtExpiresAt && jwtExpiresAt - Date.now() > 5 * 60 * 1000) {
      return jwt;
    }
    if (apiKey) {
      await get().refreshToken();
      return get().jwt;
    }
    return null;
  },
}));
