/** Format a number to fixed decimal places, stripping trailing zeros. */
export function num(value: number, decimals = 2): string {
  return parseFloat(value.toFixed(decimals)).toString();
}

/** Format USD currency. */
export function usd(value: number): string {
  return new Intl.NumberFormat("en-US", {
    style: "currency",
    currency: "USD",
    minimumFractionDigits: 0,
    maximumFractionDigits: 2,
  }).format(value);
}

/** Format percentage. */
export function pct(value: number): string {
  return `${num(value * 100, 1)}%`;
}
