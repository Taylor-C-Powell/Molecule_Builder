export function Footer() {
  return (
    <footer className="border-t border-border py-8 mt-auto">
      <div className="max-w-[1100px] mx-auto px-6 flex flex-col sm:flex-row justify-between items-center gap-4">
        <span className="text-sm text-text-muted">
          MolBuilder API &mdash; Built by Materia Foundation
        </span>
        <div className="flex gap-6">
          <a
            href="https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/docs/tutorials"
            target="_blank"
            rel="noopener noreferrer"
            className="text-sm text-text-muted hover:text-text-secondary transition-colors"
          >
            Docs
          </a>
          <a
            href="https://molbuilder-api-production.up.railway.app/docs"
            target="_blank"
            rel="noopener noreferrer"
            className="text-sm text-text-muted hover:text-text-secondary transition-colors"
          >
            API Reference
          </a>
          <a
            href="https://github.com/Taylor-C-Powell/Molecule_Builder"
            target="_blank"
            rel="noopener noreferrer"
            className="text-sm text-text-muted hover:text-text-secondary transition-colors"
          >
            GitHub
          </a>
        </div>
      </div>
    </footer>
  );
}
