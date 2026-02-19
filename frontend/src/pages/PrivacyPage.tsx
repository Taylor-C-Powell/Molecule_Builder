import { useEffect, useState } from "react";
import { Layout } from "@/components/layout/Layout";

interface Section {
  heading: string;
  content: string;
}

interface PrivacyData {
  title: string;
  effective_date: string;
  sections: Section[];
}

export default function PrivacyPage() {
  const [data, setData] = useState<PrivacyData | null>(null);

  useEffect(() => {
    fetch("/api/v1/legal/privacy")
      .then((r) => r.json())
      .then(setData)
      .catch(() => {});
  }, []);

  return (
    <Layout>
      <div className="max-w-[750px] mx-auto py-12 px-6">
        <h1 className="text-3xl font-bold mb-2">{data?.title ?? "Privacy Policy"}</h1>
        {data?.effective_date && (
          <p className="text-sm text-text-muted mb-8">
            Effective: {data.effective_date}
          </p>
        )}
        {data?.sections.map((s) => (
          <div key={s.heading} className="mb-6">
            <h2 className="text-lg font-semibold mb-2">{s.heading}</h2>
            <p className="text-sm text-text-secondary leading-relaxed">{s.content}</p>
          </div>
        ))}
        {!data && (
          <p className="text-text-muted text-sm">Loading privacy policy...</p>
        )}
      </div>
    </Layout>
  );
}
