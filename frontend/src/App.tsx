import { Routes, Route } from "react-router-dom";
import { lazy, Suspense } from "react";
import { Layout } from "./components/layout/Layout";
import { AuthGuard } from "./components/auth/AuthGuard";
import { LandingPage } from "./pages/LandingPage";

const DashboardPage = lazy(() => import("./pages/DashboardPage"));
const RetrosynthesisPage = lazy(() => import("./pages/RetrosynthesisPage"));
const ProcessPage = lazy(() => import("./pages/ProcessPage"));
const AccountPage = lazy(() => import("./pages/AccountPage"));
const RegisterPage = lazy(() => import("./pages/RegisterPage"));
const LoginPage = lazy(() => import("./pages/LoginPage"));
const LibraryPage = lazy(() => import("./pages/LibraryPage"));
const BatchPage = lazy(() => import("./pages/BatchPage"));
const TermsPage = lazy(() => import("./pages/TermsPage"));
const PrivacyPage = lazy(() => import("./pages/PrivacyPage"));

function PageLoader() {
  return (
    <div className="flex items-center justify-center min-h-[60vh]">
      <div className="h-8 w-8 animate-spin rounded-full border-2 border-accent border-t-transparent" />
    </div>
  );
}

export function App() {
  return (
    <Suspense fallback={<PageLoader />}>
      <Routes>
        <Route path="/" element={<LandingPage />} />
        <Route path="/register" element={<Layout><RegisterPage /></Layout>} />
        <Route path="/login" element={<Layout><LoginPage /></Layout>} />
        <Route path="/terms" element={<TermsPage />} />
        <Route path="/privacy" element={<PrivacyPage />} />
        <Route
          path="/dashboard"
          element={
            <AuthGuard>
              <Layout><DashboardPage /></Layout>
            </AuthGuard>
          }
        />
        <Route
          path="/retrosynthesis"
          element={
            <AuthGuard>
              <Layout><RetrosynthesisPage /></Layout>
            </AuthGuard>
          }
        />
        <Route
          path="/process"
          element={
            <AuthGuard>
              <Layout><ProcessPage /></Layout>
            </AuthGuard>
          }
        />
        <Route
          path="/library"
          element={
            <AuthGuard>
              <Layout><LibraryPage /></Layout>
            </AuthGuard>
          }
        />
        <Route
          path="/batch"
          element={
            <AuthGuard>
              <Layout><BatchPage /></Layout>
            </AuthGuard>
          }
        />
        <Route
          path="/account"
          element={
            <AuthGuard>
              <Layout><AccountPage /></Layout>
            </AuthGuard>
          }
        />
      </Routes>
    </Suspense>
  );
}
