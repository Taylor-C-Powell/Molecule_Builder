import { Link, useLocation } from "react-router-dom";
import { cn } from "@/lib/cn";
import { useAuthStore } from "@/stores/auth-store";
import { Button } from "@/components/ui/Button";
import { useState } from "react";

const NAV_LINKS = [
  { to: "/dashboard", label: "Dashboard" },
  { to: "/retrosynthesis", label: "Retrosynthesis" },
  { to: "/process", label: "Process" },
  { to: "/library", label: "Library" },
  { to: "/batch", label: "Batch" },
  { to: "/account", label: "Account" },
];

export function Header() {
  const location = useLocation();
  const { isAuthenticated, logout, email } = useAuthStore();
  const [menuOpen, setMenuOpen] = useState(false);

  return (
    <nav className="sticky top-0 z-50 bg-bg/85 backdrop-blur-xl border-b border-border">
      <div className="max-w-[1100px] mx-auto px-6 flex items-center justify-between h-[60px]">
        <Link to="/" className="flex items-center gap-2 font-bold text-lg text-text-primary hover:text-text-primary">
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="w-6 h-6">
            <circle cx="12" cy="6" r="2" /><circle cx="6" cy="18" r="2" /><circle cx="18" cy="18" r="2" />
            <line x1="12" y1="8" x2="6" y2="16" /><line x1="12" y1="8" x2="18" y2="16" /><line x1="8" y1="18" x2="16" y2="18" />
          </svg>
          MolBuilder
        </Link>

        {/* Desktop nav */}
        <div className="hidden md:flex items-center gap-6">
          {isAuthenticated ? (
            <>
              {NAV_LINKS.map((link) => (
                <Link
                  key={link.to}
                  to={link.to}
                  className={cn(
                    "text-sm font-medium transition-colors",
                    location.pathname === link.to
                      ? "text-text-primary"
                      : "text-text-secondary hover:text-text-primary",
                  )}
                >
                  {link.label}
                </Link>
              ))}
              <span className="text-xs text-text-muted">{email}</span>
              <Button variant="ghost" size="sm" onClick={logout}>
                Logout
              </Button>
            </>
          ) : (
            <>
              <Link to="/login" className="text-sm font-medium text-text-secondary hover:text-text-primary transition-colors">
                Login
              </Link>
              <Link to="/register">
                <Button size="sm">Get API Key</Button>
              </Link>
            </>
          )}
        </div>

        {/* Mobile menu button */}
        <button
          className="md:hidden text-text-secondary text-2xl bg-transparent border-none cursor-pointer"
          onClick={() => setMenuOpen(!menuOpen)}
          aria-label="Menu"
        >
          &#9776;
        </button>
      </div>

      {/* Mobile nav */}
      {menuOpen && (
        <div className="md:hidden border-b border-border bg-bg/97 p-4 flex flex-col gap-3">
          {isAuthenticated ? (
            <>
              {NAV_LINKS.map((link) => (
                <Link
                  key={link.to}
                  to={link.to}
                  onClick={() => setMenuOpen(false)}
                  className="text-sm font-medium text-text-secondary hover:text-text-primary"
                >
                  {link.label}
                </Link>
              ))}
              <Button variant="ghost" size="sm" onClick={() => { logout(); setMenuOpen(false); }}>
                Logout
              </Button>
            </>
          ) : (
            <>
              <Link to="/login" onClick={() => setMenuOpen(false)} className="text-sm text-text-secondary">Login</Link>
              <Link to="/register" onClick={() => setMenuOpen(false)}>
                <Button size="sm">Get API Key</Button>
              </Link>
            </>
          )}
        </div>
      )}
    </nav>
  );
}
