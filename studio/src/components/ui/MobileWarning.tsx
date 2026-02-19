import { useState, useEffect } from "react";

const MIN_WIDTH = 768;

export function MobileWarning() {
  const [dismissed, setDismissed] = useState(false);
  const [isMobile, setIsMobile] = useState(false);

  useEffect(() => {
    function check() {
      setIsMobile(window.innerWidth < MIN_WIDTH);
    }
    check();
    window.addEventListener("resize", check);
    return () => window.removeEventListener("resize", check);
  }, []);

  if (dismissed || !isMobile) return null;

  return (
    <div className="fixed inset-0 z-[500] flex items-center justify-center bg-bg p-6">
      <div className="text-center space-y-4 max-w-sm">
        <div className="text-lg font-bold text-text-primary">
          Desktop Recommended
        </div>
        <p className="text-sm text-text-muted">
          MolBuilder Studio is a 3D molecule editor designed for desktop
          browsers. For the best experience, please use a screen at least
          1024px wide.
        </p>
        <button
          onClick={() => setDismissed(true)}
          className="px-4 py-2 text-xs bg-accent/20 text-accent border border-accent/30 rounded-[var(--radius-sm)] hover:bg-accent/30 transition-colors cursor-pointer"
        >
          Continue Anyway
        </button>
      </div>
    </div>
  );
}
