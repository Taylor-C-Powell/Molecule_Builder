import { useEffect, useState, type ReactNode } from "react";
import { Navigate } from "react-router-dom";
import { useAuthStore } from "@/stores/auth-store";

interface AuthGuardProps {
  children: ReactNode;
}

export function AuthGuard({ children }: AuthGuardProps) {
  const isAuthenticated = useAuthStore((s) => s.isAuthenticated);
  const getValidToken = useAuthStore((s) => s.getValidToken);
  const [ready, setReady] = useState(false);

  useEffect(() => {
    if (isAuthenticated) {
      getValidToken().then(() => setReady(true));
    } else {
      setReady(true);
    }
  }, [isAuthenticated, getValidToken]);

  if (!ready) {
    return (
      <div className="h-screen flex items-center justify-center bg-bg">
        <div className="h-6 w-6 animate-spin rounded-full border-2 border-accent border-t-transparent" />
      </div>
    );
  }

  if (!isAuthenticated) {
    return <Navigate to="/login" replace />;
  }

  return <>{children}</>;
}
