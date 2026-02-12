import { Navigate } from "react-router-dom";
import { useAuthStore } from "@/stores/auth-store";
import { useEffect, useState } from "react";

interface AuthGuardProps {
  children: React.ReactNode;
}

export function AuthGuard({ children }: AuthGuardProps) {
  const { isAuthenticated, getValidToken } = useAuthStore();
  const [ready, setReady] = useState(false);

  useEffect(() => {
    if (isAuthenticated) {
      getValidToken().then(() => setReady(true));
    }
  }, [isAuthenticated, getValidToken]);

  if (!isAuthenticated) {
    return <Navigate to="/login" replace />;
  }

  if (!ready) {
    return (
      <div className="flex items-center justify-center min-h-[60vh]">
        <div className="h-8 w-8 animate-spin rounded-full border-2 border-accent border-t-transparent" />
      </div>
    );
  }

  return <>{children}</>;
}
