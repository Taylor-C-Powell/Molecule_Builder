# How I Built a Synthesis Planner You Can Pip Install

*What happens when you try to compete with C++ cheminformatics toolkits using
pure Python, stubbornness, and a lot of numpy?*

---

I started MolBuilder because I wanted to go from a SMILES string to a
manufacturing plan in one library. That tool didn't exist. RDKit is incredible
for molecular manipulation, but it doesn't plan syntheses. ASKCOS plans
syntheses, but it's a heavyweight web service you can't embed in a script.
Commercial tools like Reaxys and SciFinder cover everything, but they cost five
figures and lock you behind a browser.

I wanted something a grad student could `pip install` on a laptop and
immediately start asking: *How do I make this molecule? What reactor do I need?
What will it cost at 100 kg scale?*

So I built it. MolBuilder is a pure-Python molecular engineering toolkit --
1,288 tests, Apache 2.0, zero C++ dependencies -- that covers SMILES parsing,
3D coordinate generation, retrosynthesis planning, and process engineering. This
is the story of the three engineering decisions I found most interesting along
the way.

## The Problem with "Just Use RDKit"

RDKit is the standard for cheminformatics -- molecular descriptors, fingerprints,
substructure search, conformer generation. It's excellent at what it does. But it
doesn't do what I needed: plan a synthesis, pick a reactor, estimate costs, flag
safety hazards, and analyze scale-up. That pipeline from molecule to manufacturing
didn't exist in a single open-source package.

I also wanted something where the source code itself was a teaching tool. RDKit is
a C++ library with Python bindings -- powerful, but you can't walk a student through
the retrosynthesis scoring function by reading the source. I wanted pure Python
that a grad student could read, modify, and learn from.

MolBuilder's constraint from day one was: **numpy, scipy, matplotlib, and
nothing else.** Pure Python, readable source, minimal dependencies. This
constraint shaped every technical decision that followed.

## Decision 1: Pure-Python 3D Coordinates

The first real test of the pure-Python constraint was 3D coordinate generation.
Given a molecular graph (atoms and bonds), generate physically reasonable 3D
positions. RDKit does this with ETKDG, a mature C++ implementation. I needed to
do it in Python, fast enough to be usable.

The pipeline has five stages:

**Distance geometry embedding.** Build an N x N matrix of distance bounds --
lower and upper limits for every atom pair, derived from bond lengths, bond
angles, and topological distance. Sample random distances within those bounds,
then convert the distance matrix to a Gram matrix via double-centering and
extract 3D coordinates from the top three eigenvalues. This is textbook distance
geometry, but the implementation details matter: you need `numpy.linalg.eigh`
for the eigendecomposition and careful handling of negative eigenvalues (which
happen when your sampled distances are geometrically inconsistent).

**Violation minimization.** The initial embedding usually has some distance
bounds violations. L-BFGS-B from scipy minimizes a squared-violation penalty
function. This gets you coordinates that respect the distance bounds but look
nothing like real molecules yet -- bond lengths are right, but rings are mangled
and angles are off.

**Ring template correction.** This was the hardest part. Rings need to look like
rings -- planar pentagons, cyclohexane chairs, aromatic hexagons. I detect rings
up to size 8 with DFS, generate ideal templates (regular polygons for planar
rings, puckered geometries for saturated rings), then align each template to the
embedded positions using Kabsch SVD alignment. Fused rings are trickier: they
share edges, so you anchor on the shared atoms and rotate the second ring into
place.

**Force field optimization.** A simple forcefield with bond stretching, angle
bending, torsion rotation, and Lennard-Jones nonbonded terms. Equilibrium values
come from hybridization: SP3 gets 109.47 degrees, SP2 gets 120, SP gets 180.
SP2 centers get an extra out-of-plane penalty to keep them flat. Optimized with
L-BFGS-B again.

**Stereo enforcement.** Check the signed volume at each chiral center and swap
atoms if the chirality is wrong.

The result is not as good as ETKDG. The conformers are reasonable -- bonds and
angles look right, rings are recognizable, chirality is correct -- but the
torsional sampling is limited and the forcefield is basic. For visualization and
rough geometry it works well. For publication-quality conformer generation, use
RDKit. MolBuilder supports RDKit as an optional backend (`pip install
molbuilder[rdkit]`) for exactly this reason.

**What I learned:** Competing with C++ on numerical quality is a losing game.
Competing on *accessibility* is a winning one. Most users need "good enough"
coordinates to see what their molecule looks like, not Boltzmann-weighted
conformer ensembles. The pure-Python implementation serves that use case without
any installation friction.

## Decision 2: Scoring Retrosynthesis Like a Chemist

MolBuilder's retrosynthesis engine uses beam search over 185 reaction templates
across 14 categories (substitution, coupling, carbonyl, pericyclic, etc.). Given
a target molecule, it proposes disconnections, scores them, and keeps the top k
at each level. The tree terminates when all leaves are purchasable starting
materials (from a knowledge base of 270 compounds) or the depth limit is
reached.

The interesting part is the scoring function. Each disconnection gets a score
from 0 to 100 based on five factors:

1. **Yield expectation** (0-25 pts): The midpoint of the template's typical
   yield range. A Suzuki coupling (70-95%) scores higher than a radical
   bromination (40-60%).

2. **Precursor availability** (0-30 pts): What fraction of the precursors are in
   the purchasable materials database. A disconnection that produces two
   purchasable fragments scores 30. One that produces an exotic intermediate
   scores lower.

3. **Complexity reduction** (0-20 pts): Does the disconnection actually simplify
   the problem? Measured by heavy atom count -- if the largest precursor has
   fewer heavy atoms than the target, you're making progress. There's a
   mass-balance penalty too: if total precursor atoms are less than 40% of the
   target, the disconnection is probably losing important structure.

4. **Strategic bond preference** (0-15 pts): C-C bond forming reactions score
   highest. Coupling reactions (Suzuki, Heck, Sonogashira, Stille, Negishi) get
   the maximum bonus. This encodes the synthetic chemistry heuristic that you
   should plan your C-C bond formations first and save functional group
   manipulations for later.

5. **Category bonus** (0-10 pts): A tiebreaker that favors coupling and carbonyl
   chemistry over protection/deprotection steps.

These weights encode decades of synthetic chemistry intuition. The strategic bond
preference, for instance, reflects Corey's retrosynthetic analysis principles:
disconnect at strategic bonds (usually C-C) that simplify the molecule into
roughly equal-sized fragments. The availability weighting reflects the practical
reality that a beautiful retrosynthetic plan is useless if you can't buy the
starting materials.

Is this as good as a neural retrosynthesis model trained on millions of
reactions? No. But it's transparent -- you can read the scoring function and
understand exactly why a route was chosen. When a student asks "why did it pick
Suzuki coupling here instead of Grignard addition?", the answer is in the
weights, not in a black box.

## Decision 3: Making Template Matching Chemistry-Aware

The initial retrosynthesis engine matched templates purely on functional group
patterns. If your molecule had a ketone and an amine, it would suggest reductive
amination. Simple and fast, but chemistry-blind: it couldn't tell the difference
between an unhindered primary amine (fast reaction, high yield) and a
sterically-blocked tertiary amine next to a neopentyl group (good luck).

The condition prediction module adds substrate awareness through three
classification layers:

**Steric classification.** For each functional group center, count the branching
degree of the atom itself plus its alpha carbons. Unhindered (< 4 total branch
points), moderately hindered (4-5), or hindered (>= 6). This matters enormously
for reaction selection: SN2 displacement on a neopentyl halide is essentially
impossible, and the scoring function now knows that. Hindered substrates get
zero points for SN2-type templates.

**Electronic character.** Count electron-donating groups (alcohols, ethers,
amines, alkenes) and electron-withdrawing groups (nitriles, nitros, carbonyls,
sulfonamides). Aromatic rings contribute +0.5 EDG each. The balance determines
if the substrate is electron-rich, neutral, or electron-poor. Electrophilic
aromatic substitution on a nitrobenzene? Low score. On an anisole? High score.

**Sensitive group identification.** Flag functional groups that might not survive
the reaction conditions -- aldehydes, epoxides, free thiols, boronic acids. The
scoring function applies a compatibility penalty (up to -40 points) when the
proposed conditions would destroy a sensitive group elsewhere in the molecule.

This isn't a quantum mechanical calculation. It's pattern matching with chemical
knowledge encoded in the patterns. But it catches the most common mistakes
students make when proposing synthetic routes: ignoring sterics, ignoring
electronics, and ignoring functional group compatibility.

## The Numbers

As of v1.2.0:

- 185 reaction templates across 14 categories
- 24 functional group detectors
- 270 purchasable starting materials
- 1,288 tests across 24 test files
- Pure Python: numpy, scipy, matplotlib
- Optional RDKit backend for higher-quality 3D
- 9 tutorials (3 Jupyter notebooks + 6 written guides)
- SaaS API with free tier (100 requests/day)
- Apache 2.0 license

## What I Got Wrong

Honesty section.

**The template library is small.** 185 templates sounds like a lot until you
compare it to Reaxys (70M+ reactions) or even ASKCOS (~160k templates). For
common organic chemistry -- undergrad through first-year grad -- the coverage is
decent. For niche transformations, natural product synthesis, or anything
involving exotic metals, you'll hit gaps quickly.

**Rule-based retrosynthesis has a ceiling.** Neural models trained on large
reaction databases produce better routes for complex targets. MolBuilder's beam
search with hand-tuned weights is interpretable and predictable, but it won't
surprise you with a creative disconnection the way a learned model might.

**The forcefield is basic.** No explicit electrostatics, no solvation, limited
torsion parameterization. Fine for visualization and rough geometry. Not suitable
for energy calculations, docking, or anything quantitative.

**The purchasable materials list is curated, not live.** 270 compounds is a
reasonable starting set, but real availability checking requires querying vendor
APIs (Sigma-Aldrich, TCI, Enamine). That's on the roadmap.

These are real limitations, not "future work" euphemisms. If your use case bumps
into them, MolBuilder is the wrong tool today. I'd rather you know that upfront.

## What's Next

The SaaS API is live with a free academic tier. The immediate priorities are
expanding the template library (community contributions welcome), adding vendor
API integration for real-time availability, and building a proper web-based
route visualization.

Longer term, I want to explore hybrid approaches -- using the rule-based engine
for route scoring and feasibility assessment while delegating creative
disconnection proposals to a learned model. The scoring function and condition
prediction are MolBuilder's real strengths; the template matching is the
bottleneck.

If you're teaching computational chemistry, process chemistry, or Python for
scientists, I'd especially love feedback on what's useful and what's missing. The
tutorials were designed for classroom use, and I'm actively looking for
professors willing to pilot them.

---

**Try it:**

```
pip install molbuilder
```

**Links:**
- GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
- PyPI: https://pypi.org/project/molbuilder/
- Tutorials: https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/tutorials
- API: https://molbuilder-api.up.railway.app
