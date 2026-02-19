import { useState } from "react";
import { useNavigate, Link } from "react-router-dom";
import { useAuthStore } from "@/stores/auth-store";
import { Button } from "@/components/ui/Button";
import { Input } from "@/components/ui/Input";
import { Alert } from "@/components/ui/Alert";

export function LoginForm() {
  const [apiKey, setApiKey] = useState("");
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const navigate = useNavigate();
  const login = useAuthStore((s) => s.login);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    if (!apiKey.trim()) return;
    setLoading(true);
    setError(null);

    try {
      await login(apiKey.trim());
      navigate("/");
    } catch (err) {
      setError(err instanceof Error ? err.message : "Login failed");
    } finally {
      setLoading(false);
    }
  }

  return (
    <form onSubmit={handleSubmit} className="space-y-4 max-w-sm mx-auto">
      <h1 className="text-xl font-bold text-text-primary">Sign in to Studio</h1>
      <p className="text-sm text-text-secondary">
        Enter your MolBuilder API key to access the studio.
      </p>

      {error && <Alert variant="error">{error}</Alert>}

      <Input
        label="API Key"
        type="password"
        placeholder="mb_..."
        value={apiKey}
        onChange={(e) => setApiKey(e.target.value)}
      />

      <Button type="submit" disabled={loading || !apiKey.trim()} className="w-full">
        {loading ? "Signing in..." : "Sign In"}
      </Button>

      <p className="text-xs text-text-muted text-center">
        Don't have an API key?{" "}
        <Link to="/register" className="text-accent hover:underline">
          Register
        </Link>
      </p>
    </form>
  );
}
