import { Link } from "react-router-dom";
import { Button } from "@/components/ui/Button";
import { Footer } from "@/components/layout/Footer";
import { useAuthStore } from "@/stores/auth-store";
import { useState } from "react";

const FEATURES = [
  {
    title: "SMILES Parsing",
    desc: "Parse any SMILES string into a full molecular graph with atom counts, bond topology, and canonical normalization.",
    endpoint: "POST /api/v1/molecule/from-smiles",
    icon: (
      <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
        <circle cx="12" cy="12" r="3" /><circle cx="5" cy="6" r="2" /><circle cx="19" cy="6" r="2" />
        <circle cx="5" cy="18" r="2" /><circle cx="19" cy="18" r="2" />
        <line x1="7" y1="7" x2="10" y2="10" /><line x1="14" y1="10" x2="17" y2="7" />
        <line x1="7" y1="17" x2="10" y2="14" /><line x1="14" y1="14" x2="17" y2="17" />
      </svg>
    ),
  },
  {
    title: "Retrosynthetic Analysis",
    desc: "Route planning with purchasable starting materials, yield predictions, and named reaction identification.",
    endpoint: "POST /api/v1/retrosynthesis/plan",
    icon: (
      <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
        <polyline points="22 12 18 12 15 21 9 3 6 12 2 12" />
      </svg>
    ),
  },
  {
    title: "Process Engineering",
    desc: "Safety analysis, cost breakdown, reactor sizing, purification steps, and scale-up recommendations.",
    endpoint: "POST /api/v1/process/evaluate",
    icon: (
      <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
        <rect x="2" y="7" width="20" height="14" rx="2" />
        <path d="M16 7V4a2 2 0 0 0-2-2h-4a2 2 0 0 0-2 2v3" />
        <line x1="12" y1="11" x2="12" y2="17" /><line x1="9" y1="14" x2="15" y2="14" />
      </svg>
    ),
  },
  {
    title: "3D Coordinates",
    desc: "Optimized 3D atomic positions with hybridization states, formal charges, and bond geometry.",
    endpoint: "GET /api/v1/molecule/{id}/3d",
    icon: (
      <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
        <path d="M12 3v18" /><path d="M3 12h18" /><path d="M3.5 5.5l17 13" />
      </svg>
    ),
  },
  {
    title: "Molecular Properties",
    desc: "Formula, molecular weight, Lipinski analysis, functional group detection, logP, TPSA, and more.",
    endpoint: "GET /api/v1/molecule/{id}/properties",
    icon: (
      <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
        <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z" />
        <polyline points="14 2 14 8 20 8" /><line x1="9" y1="15" x2="15" y2="15" />
      </svg>
    ),
  },
  {
    title: "Element Lookup",
    desc: "Instant periodic table data: atomic number, weight, symbol, and element name for any element.",
    endpoint: "GET /api/v1/elements/{symbol}",
    icon: (
      <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
        <rect x="3" y="3" width="7" height="7" /><rect x="14" y="3" width="7" height="7" />
        <rect x="3" y="14" width="7" height="7" /><rect x="14" y="14" width="7" height="7" />
      </svg>
    ),
  },
];

export function LandingPage() {
  const isAuthenticated = useAuthStore((s) => s.isAuthenticated);
  const [menuOpen, setMenuOpen] = useState(false);

  return (
    <div className="min-h-screen bg-bg text-text-primary">
      {/* Nav */}
      <nav className="sticky top-0 z-50 bg-bg/85 backdrop-blur-xl border-b border-border">
        <div className="max-w-[1100px] mx-auto px-6 flex items-center justify-between h-[60px]">
          <Link to="/" className="flex items-center gap-2 font-bold text-lg text-text-primary">
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="w-6 h-6">
              <circle cx="12" cy="6" r="2" /><circle cx="6" cy="18" r="2" /><circle cx="18" cy="18" r="2" />
              <line x1="12" y1="8" x2="6" y2="16" /><line x1="12" y1="8" x2="18" y2="16" /><line x1="8" y1="18" x2="16" y2="18" />
            </svg>
            MolBuilder
          </Link>
          <div className="hidden md:flex items-center gap-8">
            <a href="#features" className="text-sm font-medium text-text-secondary hover:text-text-primary transition-colors">Features</a>
            <a href="#pricing" className="text-sm font-medium text-text-secondary hover:text-text-primary transition-colors">Pricing</a>
            <a href="https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/docs/tutorials" target="_blank" rel="noopener noreferrer" className="text-sm font-medium text-text-secondary hover:text-text-primary transition-colors">Docs</a>
            {isAuthenticated ? (
              <Link to="/dashboard"><Button size="sm">Dashboard</Button></Link>
            ) : (
              <Link to="/register"><Button size="sm">Get API Key</Button></Link>
            )}
          </div>
          <button className="md:hidden text-text-secondary text-2xl bg-transparent border-none cursor-pointer" onClick={() => setMenuOpen(!menuOpen)} aria-label="Menu">&#9776;</button>
        </div>
        {menuOpen && (
          <div className="md:hidden border-b border-border bg-bg/97 p-4 flex flex-col gap-3">
            <a href="#features" onClick={() => setMenuOpen(false)} className="text-sm text-text-secondary">Features</a>
            <a href="#pricing" onClick={() => setMenuOpen(false)} className="text-sm text-text-secondary">Pricing</a>
            {isAuthenticated ? (
              <Link to="/dashboard" onClick={() => setMenuOpen(false)}><Button size="sm">Dashboard</Button></Link>
            ) : (
              <Link to="/register" onClick={() => setMenuOpen(false)}><Button size="sm">Get API Key</Button></Link>
            )}
          </div>
        )}
      </nav>

      {/* Hero */}
      <section className="relative py-24 md:py-32 text-center overflow-hidden">
        <div className="absolute top-[-40%] left-1/2 -translate-x-1/2 w-[800px] h-[600px] bg-[radial-gradient(ellipse,var(--color-accent-glow)_0%,transparent_70%)] pointer-events-none" />
        <div className="max-w-[1100px] mx-auto px-6 relative">
          <h1 className="text-4xl md:text-6xl font-extrabold tracking-tight leading-tight mb-5">
            Molecular Engineering Platform
          </h1>
          <p className="text-lg md:text-xl text-text-secondary max-w-[600px] mx-auto mb-10">
            Parse SMILES strings, visualize molecules in 3D, plan retrosynthesis routes, and evaluate manufacturing processes.
          </p>
          <div className="flex gap-4 justify-center flex-wrap">
            {isAuthenticated ? (
              <Link to="/dashboard"><Button size="lg">Open Dashboard</Button></Link>
            ) : (
              <Link to="/register"><Button size="lg">Get API Key &mdash; Free</Button></Link>
            )}
            <a href="https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/docs/tutorials" target="_blank" rel="noopener noreferrer">
              <Button variant="secondary" size="lg">View Documentation</Button>
            </a>
          </div>
        </div>
      </section>

      {/* Features */}
      <section className="border-t border-border py-20" id="features">
        <div className="max-w-[1100px] mx-auto px-6">
          <h2 className="text-3xl font-bold tracking-tight text-center mb-2">API Capabilities</h2>
          <p className="text-center text-text-secondary mb-12">Everything you need for computational chemistry, accessible over HTTP.</p>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-5">
            {FEATURES.map((f) => (
              <div key={f.title} className="p-6 bg-bg-card border border-border rounded-[var(--radius-md)] hover:border-border-hover transition-colors">
                <div className="w-10 h-10 bg-accent-glow border border-accent/20 rounded-[var(--radius-md)] flex items-center justify-center text-accent mb-4">
                  {f.icon}
                </div>
                <h3 className="font-semibold mb-2">{f.title}</h3>
                <p className="text-sm text-text-secondary mb-3">{f.desc}</p>
                <code className="text-xs text-text-muted font-mono">{f.endpoint}</code>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* Pricing */}
      <section className="border-t border-border py-20" id="pricing">
        <div className="max-w-[1100px] mx-auto px-6">
          <h2 className="text-3xl font-bold tracking-tight text-center mb-2">Pricing</h2>
          <p className="text-center text-text-secondary mb-12">Start free. Upgrade when you need more throughput.</p>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 max-w-[750px] mx-auto">
            {/* Free */}
            <div className="p-8 bg-bg-card border border-border rounded-[var(--radius-md)] text-center">
              <h3 className="text-xl font-semibold mb-3">Free</h3>
              <div className="text-5xl font-extrabold tracking-tight mb-1">$0<span className="text-base font-normal text-text-muted">/month</span></div>
              <p className="text-sm text-text-muted mb-6">No credit card required</p>
              <ul className="text-left space-y-2 mb-8">
                {["10 requests per minute", "5 expensive operations per hour", "All endpoints included", "Community support"].map((item) => (
                  <li key={item} className="flex items-center gap-2 text-sm text-text-secondary">
                    <svg className="w-4 h-4 text-green flex-shrink-0" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5"><path strokeLinecap="round" strokeLinejoin="round" d="M5 13l4 4L19 7" /></svg>
                    {item}
                  </li>
                ))}
              </ul>
              <Link to="/register" className="block">
                <Button variant="secondary" className="w-full">Get Started</Button>
              </Link>
            </div>
            {/* Pro */}
            <div className="relative p-8 bg-bg-card border border-accent rounded-[var(--radius-md)] text-center shadow-[0_0_50px_var(--color-accent-glow)]">
              <span className="absolute -top-3 left-1/2 -translate-x-1/2 bg-accent text-white text-[0.7rem] font-bold uppercase tracking-wider px-3 py-1 rounded-full">Popular</span>
              <h3 className="text-xl font-semibold mb-3">Pro</h3>
              <div className="text-5xl font-extrabold tracking-tight mb-1">$49<span className="text-base font-normal text-text-muted">/month</span></div>
              <p className="text-sm text-text-muted mb-6">or $490/year (save 17%)</p>
              <ul className="text-left space-y-2 mb-8">
                {["60 requests per minute", "30 expensive operations per hour", "All endpoints included", "Priority support"].map((item) => (
                  <li key={item} className="flex items-center gap-2 text-sm text-text-secondary">
                    <svg className="w-4 h-4 text-green flex-shrink-0" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5"><path strokeLinecap="round" strokeLinejoin="round" d="M5 13l4 4L19 7" /></svg>
                    {item}
                  </li>
                ))}
              </ul>
              <Link to="/register" className="block">
                <Button className="w-full">Upgrade to Pro</Button>
              </Link>
            </div>
          </div>
        </div>
      </section>

      <Footer />
    </div>
  );
}
