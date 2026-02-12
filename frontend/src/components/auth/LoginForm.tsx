import { useState } from "react";
import { useNavigate, Link } from "react-router-dom";
import { useAuthStore } from "@/stores/auth-store";
import { Button } from "@/components/ui/Button";
import { Input } from "@/components/ui/Input";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Alert } from "@/components/ui/Alert";

export function LoginForm() {
  const [apiKey, setApiKey] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const { setApiKey: storeKey, login } = useAuthStore();
  const navigate = useNavigate();

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    setLoading(true);
    setError(null);
    try {
      // Verify the key works by exchanging for JWT
      await login(apiKey);
      storeKey(apiKey, "", "free");
      navigate("/dashboard");
    } catch (err) {
      setError(err instanceof Error ? err.message : "Invalid API key");
    } finally {
      setLoading(false);
    }
  }

  return (
    <Card className="max-w-md mx-auto mt-12">
      <CardHeader>
        <CardTitle>Sign In</CardTitle>
        <p className="text-sm text-text-secondary mt-1">
          Enter your API key to access the dashboard.
        </p>
      </CardHeader>
      <form onSubmit={handleSubmit} className="flex flex-col gap-4">
        {error && <Alert variant="error">{error}</Alert>}
        <Input
          label="API Key"
          type="password"
          value={apiKey}
          onChange={(e) => setApiKey(e.target.value)}
          placeholder="mb_..."
          required
          autoFocus
        />
        <Button type="submit" disabled={loading || !apiKey} className="w-full">
          {loading ? "Verifying..." : "Sign In"}
        </Button>
        <p className="text-sm text-text-muted text-center">
          Don't have a key?{" "}
          <Link to="/register" className="text-accent hover:text-accent-hover">
            Register for free
          </Link>
        </p>
      </form>
    </Card>
  );
}
