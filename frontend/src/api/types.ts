// ---- Auth ----
export interface RegisterRequest {
  email: string;
}

export interface RegisterResponse {
  api_key: string;
  email: string;
  tier: string;
  role: string;
  message: string;
}

// Note: /auth/token accepts API key in X-API-Key header, not body
export interface TokenRequest {
  api_key: string;
}

export interface TokenResponse {
  access_token: string;
  token_type: string;
  expires_in: number;
  email: string;
  tier: string;
}

// ---- Molecule ----
export interface ParseSmilesRequest {
  smiles: string;
  name?: string;
}

export interface MoleculeResponse {
  id: string;
  name: string;
  smiles: string;
  num_atoms: number;
  num_bonds: number;
}

export interface AtomResponse {
  index: number;
  symbol: string;
  position: [number, number, number];
  hybridization: string | null;
  formal_charge: number;
}

export interface BondResponse {
  atom_i: number;
  atom_j: number;
  order: number;
  rotatable: boolean;
}

export interface Molecule3DResponse {
  id: string;
  atoms: AtomResponse[];
  bonds: BondResponse[];
}

export interface MoleculePropertiesResponse {
  id: string;
  smiles: string;
  formula: string;
  molecular_weight: number;
  num_atoms: number;
  num_bonds: number;
  functional_groups: string[];
  logp: number | null;
  hbd: number | null;
  hba: number | null;
  rotatable_bonds: number | null;
  tpsa: number | null;
  heavy_atom_count: number | null;
  lipinski_violations: number | null;
  lipinski_pass: boolean | null;
  sa_score: number | null;
}

// ---- ADMET ----
export interface ADMETResponse {
  id: string;
  smiles: string;
  oral_bioavailability: string;
  intestinal_absorption: string;
  caco2_permeability: string;
  pgp_substrate: boolean;
  bbb_penetrant: boolean;
  plasma_protein_binding: string;
  vd_class: string;
  cyp_inhibition: Record<string, boolean>;
  metabolic_stability: string;
  renal_clearance: string;
  half_life_class: string;
  herg_risk: string;
  ames_mutagenicity: boolean;
  hepatotoxicity_risk: string;
  structural_alerts: string[];
  overall_score: number;
  warnings: string[];
  flags: string[];
}

// ---- Solubility ----
export interface SolubilityResponse {
  id: string;
  smiles: string;
  log_s_esol: number;
  log_s_gse: number;
  solubility_mg_ml: number;
  solubility_class: string;
  estimated_melting_point_c: number;
  crystallization_risk: string;
  polymorph_risk: string;
}

// ---- Retrosynthesis ----
export interface RetroRequest {
  smiles: string;
  max_depth?: number;
  beam_width?: number;
}

export interface PrecursorResponse {
  smiles: string;
  name: string;
  cost_per_kg: number;
}

export interface DisconnectionResponse {
  reaction_name: string;
  named_reaction: string | null;
  category: string;
  score: number;
  precursors: PrecursorResponse[];
}

export interface RetroNodeResponse {
  smiles: string;
  is_purchasable: boolean;
  functional_groups: string[];
  best_disconnection: DisconnectionResponse | null;
  children: RetroNodeResponse[];
}

export interface SynthesisStepResponse {
  step_number: number;
  reaction_name: string;
  named_reaction: string | null;
  category: string;
  precursor_smiles: string[];
  product_smiles: string;
  product_name: string;
  conditions: string;
  expected_yield: number;
  notes: string;
}

export interface SynthesisRouteResponse {
  target_smiles: string;
  target_name: string;
  total_steps: number;
  overall_yield: number;
  longest_linear_sequence: number;
  starting_materials: PrecursorResponse[];
  steps: SynthesisStepResponse[];
}

export interface RetroResponse {
  tree: RetroNodeResponse;
  routes_found: number;
  max_depth: number;
  beam_width: number;
  best_route: SynthesisRouteResponse | null;
}

// ---- Process Engineering ----
export interface ProcessEvaluateRequest {
  smiles: string;
  scale_kg?: number;
  max_depth?: number;
  beam_width?: number;
}

export interface ReactorResponse {
  reactor_type: string;
  volume_L: number;
  temperature_C: number;
  pressure_atm: number;
  residence_time_min: number;
  mixing_type: string;
  heat_transfer: string;
  material: string;
  estimated_cost_usd: number;
  notes: string;
}

export interface ConditionsResponse {
  temperature_C: number;
  pressure_atm: number;
  solvent: string;
  concentration_M: number;
  addition_rate: string;
  reaction_time_hours: number;
  atmosphere: string;
  workup_procedure: string;
  notes: string;
}

export interface PurificationStepResponse {
  method: string;
  description: string;
  estimated_recovery: number;
  estimated_purity: number;
  scale_appropriate: boolean;
  notes: string;
}

export interface HazardInfoResponse {
  reagent_name: string;
  ghs_hazards: string[];
  ghs_pictograms: string[];
  hazard_descriptions: string[];
  pictogram_descriptions: string[];
}

export interface SafetyAssessmentResponse {
  step_number: number;
  step_name: string;
  risk_level: string;
  ppe_required: string[];
  engineering_controls: string[];
  emergency_procedures: string[];
  incompatible_materials: string[];
  waste_classification: string;
  hazards: HazardInfoResponse[];
}

export interface CostBreakdownResponse {
  raw_materials_usd: number;
  labor_usd: number;
  equipment_usd: number;
  energy_usd: number;
  waste_disposal_usd: number;
  overhead_usd: number;
}

export interface CostEstimateResponse {
  total_usd: number;
  per_kg_usd: number;
  scale_kg: number;
  breakdown: CostBreakdownResponse;
  notes: string[];
}

export interface ScaleUpResponse {
  target_annual_kg: number;
  recommended_mode: string;
  batch_size_kg: number | null;
  batches_per_year: number | null;
  cycle_time_hours: number;
  annual_capacity_kg: number;
  capital_cost_usd: number;
  operating_cost_annual_usd: number;
  scale_up_risks: string[];
  recommendations: string[];
}

export interface StepProcessDetail {
  step_number: number;
  reaction_name: string;
  reactor: ReactorResponse;
  conditions: ConditionsResponse;
  purification: PurificationStepResponse[];
}

export interface ProcessEvaluateResponse {
  smiles: string;
  scale_kg: number;
  route_found: boolean;
  total_steps: number;
  overall_yield: number;
  step_details: StepProcessDetail[];
  safety: SafetyAssessmentResponse[];
  cost: CostEstimateResponse | null;
  scale_up: ScaleUpResponse | null;
}

// ---- Billing ----
export interface CheckoutRequest {
  plan: "pro_monthly" | "pro_yearly" | "team_monthly" | "team_yearly";
  success_url: string;
  cancel_url: string;
}

export interface CheckoutResponse {
  checkout_url: string;
}

export interface PortalResponse {
  portal_url: string;
}

export interface BillingStatusResponse {
  email: string;
  tier: string;
  subscription_status: string;
  has_billing: boolean;
  stripe_customer_id: string | null;
  stripe_subscription_id: string | null;
}

// ---- Library ----
export interface LibrarySaveRequest {
  smiles: string;
  name?: string | null;
  tags?: string[];
  notes?: string | null;
}

export interface LibraryUpdateRequest {
  name?: string | null;
  tags?: string[] | null;
  notes?: string | null;
}

export interface LibraryMoleculeResponse {
  id: number;
  smiles: string;
  name: string | null;
  tags: string[];
  notes: string | null;
  properties: Record<string, unknown>;
  created_at: string;
  updated_at: string;
}

export interface LibraryListResponse {
  molecules: LibraryMoleculeResponse[];
  total: number;
  page: number;
  per_page: number;
}

export interface LibraryImportRequest {
  smiles_list: string[];
  tag?: string | null;
}

export interface LibraryImportResponse {
  saved: number;
  duplicates: number;
  errors: string[];
}

// ---- Batch ----
export interface BatchSubmitRequest {
  smiles_list: string[];
  job_type: "properties" | "retrosynthesis" | "conditions" | "evaluate";
  params?: Record<string, unknown> | null;
}

export interface BatchSubmitResponse {
  job_id: string;
  status: string;
  created_at: string;
}

export interface BatchStatusResponse {
  job_id: string;
  status: string;
  job_type: string;
  progress_pct: number;
  result: Record<string, unknown> | null;
  error: string | null;
  created_at: string;
  updated_at: string;
}

export interface BatchJobSummary {
  job_id: string;
  status: string;
  job_type: string;
  progress_pct: number;
  created_at: string;
  updated_at: string;
}

export interface BatchListResponse {
  jobs: BatchJobSummary[];
  total: number;
  page: number;
  per_page: number;
}

// ---- File I/O ----
export interface FileImportResponse {
  molecules: MoleculeResponse[];
  format: string;
  count: number;
}

// ---- Common ----
export interface ApiError {
  detail: string;
}
