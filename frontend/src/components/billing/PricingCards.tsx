import { useState } from "react";
import { Button } from "@/components/ui/Button";

interface PricingCardsProps {
  onCheckout: (plan: "pro_monthly" | "pro_yearly") => void;
  currentTier?: string;
}

const CHECK = (
  <svg className="w-4 h-4 text-green flex-shrink-0" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5">
    <path strokeLinecap="round" strokeLinejoin="round" d="M5 13l4 4L19 7" />
  </svg>
);

export function PricingCards({ onCheckout, currentTier }: PricingCardsProps) {
  const isPro = currentTier === "pro";
  const [yearly, setYearly] = useState(false);

  return (
    <div className="flex flex-col items-center gap-6">
      {/* Monthly / Yearly toggle */}
      <div className="inline-flex items-center rounded-full bg-bg-card border border-border p-1 text-sm">
        <button
          type="button"
          className={`px-4 py-1.5 rounded-full font-medium transition-colors ${
            !yearly ? "bg-accent text-white" : "text-text-secondary hover:text-text-primary"
          }`}
          onClick={() => setYearly(false)}
        >
          Monthly
        </button>
        <button
          type="button"
          className={`px-4 py-1.5 rounded-full font-medium transition-colors ${
            yearly ? "bg-accent text-white" : "text-text-secondary hover:text-text-primary"
          }`}
          onClick={() => setYearly(true)}
        >
          Yearly
        </button>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-6 max-w-[750px]">
        {/* Free */}
        <div className="p-8 bg-bg-card border border-border rounded-[var(--radius-md)] text-center">
          <h3 className="text-xl font-semibold mb-3">Free</h3>
          <div className="text-5xl font-extrabold tracking-tight mb-1">
            $0<span className="text-base font-normal text-text-muted">/month</span>
          </div>
          <p className="text-sm text-text-muted mb-6">No credit card required</p>
          <ul className="text-left space-y-2 mb-8">
            {["10 requests per minute", "5 expensive operations per hour", "All endpoints included", "Community support"].map((item) => (
              <li key={item} className="flex items-center gap-2 text-sm text-text-secondary">{CHECK} {item}</li>
            ))}
          </ul>
          <Button variant="secondary" className="w-full" disabled={!isPro}>
            {currentTier === "free" ? "Current Plan" : "Free Tier"}
          </Button>
        </div>

        {/* Pro */}
        <div className="relative p-8 bg-bg-card border border-accent rounded-[var(--radius-md)] text-center shadow-[0_0_50px_var(--color-accent-glow)]">
          <span className="absolute -top-3 left-1/2 -translate-x-1/2 bg-accent text-white text-[0.7rem] font-bold uppercase tracking-wider px-3 py-1 rounded-full">
            Popular
          </span>
          <h3 className="text-xl font-semibold mb-3">Pro</h3>
          <div className="text-5xl font-extrabold tracking-tight mb-1">
            {yearly ? "$41" : "$49"}
            <span className="text-base font-normal text-text-muted">/month</span>
          </div>
          <p className="text-sm text-text-muted mb-6">
            {yearly ? "$490/year (save 17%)" : "or $490/year (save 17%)"}
          </p>
          <ul className="text-left space-y-2 mb-8">
            {["60 requests per minute", "30 expensive operations per hour", "All endpoints included", "Priority support"].map((item) => (
              <li key={item} className="flex items-center gap-2 text-sm text-text-secondary">{CHECK} {item}</li>
            ))}
          </ul>
          <Button
            className="w-full"
            disabled={isPro}
            onClick={() => onCheckout(yearly ? "pro_yearly" : "pro_monthly")}
          >
            {isPro ? "Current Plan" : "Upgrade to Pro"}
          </Button>
        </div>
      </div>
    </div>
  );
}
