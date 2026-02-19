import { Component, type ReactNode } from "react";

interface Props {
  children: ReactNode;
  fallback?: ReactNode;
}

interface State {
  error: Error | null;
}

export class ErrorBoundary extends Component<Props, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State {
    return { error };
  }

  handleRetry = () => {
    this.setState({ error: null });
  };

  render() {
    if (this.state.error) {
      if (this.props.fallback) return this.props.fallback;

      return (
        <div className="flex items-center justify-center h-full bg-bg p-8">
          <div className="text-center space-y-3 max-w-md">
            <div className="text-sm font-semibold text-red-400">
              Something went wrong
            </div>
            <div className="text-xs text-text-muted font-mono break-all">
              {this.state.error.message}
            </div>
            <button
              onClick={this.handleRetry}
              className="px-4 py-1.5 text-xs bg-accent/20 text-accent border border-accent/30 rounded-[var(--radius-sm)] hover:bg-accent/30 transition-colors cursor-pointer"
            >
              Retry
            </button>
          </div>
        </div>
      );
    }

    return this.props.children;
  }
}
