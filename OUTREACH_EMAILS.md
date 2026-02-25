# Professor Outreach Emails

Ready-to-send email templates for MolBuilder academic outreach. Three variants
targeting different professor profiles. Personalize each before sending.

**Key links to include:**
- GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
- PyPI: https://pypi.org/project/molbuilder/
- Tutorials: https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/tutorials
- API docs: https://molbuilder-api.up.railway.app

---

## Template A -- Retrosynthesis / Process Chemistry Professors

**Best for:** Faculty teaching synthesis planning, process chemistry, or chemical
engineering courses. Highest product-curriculum fit.

**Targets:**

| Name | University | Why they fit |
|------|-----------|-------------|
| Vincent Scalfani | U. Alabama | Cheminformatics librarian, builds open chemistry data tools and workshops |
| Tim Cernak | U. Michigan | Runs Cernak Lab -- reaction informatics, synthesis automation, ML for chemistry |
| Robert Belford | UA Little Rock | Champions open-source chemistry education, OLCC (Online Chemistry Learning Community) |

---

**Subject:** Open-source retrosynthesis + process engineering toolkit for [COURSE NAME]

**Body:**

Hi Professor [LAST NAME],

I came across [PERSONALIZE: their specific course, workshop, paper, or tool -- e.g.
"your cheminformatics workshop series at Alabama" / "the Cernak Lab's work on
reaction informatics" / "your OLCC initiative for open chemistry education"] and
thought MolBuilder might be useful for your students.

MolBuilder is an open-source Python toolkit that covers the full pipeline from
molecular structure through retrosynthesis to process engineering -- something I
couldn't find in a single open package when I needed it. The highlights:

- **185 reaction templates** across 14 categories with beam-search retrosynthesis
- **Substrate-aware condition prediction** (steric/electronic classification)
- **Process engineering pipeline**: reactor selection, safety assessment, cost
  estimation, scale-up analysis
- Pure Python -- `pip install molbuilder`, no RDKit required
- 1,280+ tests, Apache 2.0 licensed

I've included Jupyter tutorials that walk through SMILES-to-synthesis and
process engineering workflows:
https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/tutorials

Would you be open to taking a look? I'd value any feedback, and if it fits your
curriculum, I'd be glad to help adapt examples for [PERSONALIZE: their specific
course context].

Best,
Taylor Powell

GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
PyPI: https://pypi.org/project/molbuilder/

---

### Personalization Notes for Template A

- **Scalfani:** Reference his cheminformatics workshops and open data advocacy at
  Alabama Libraries. He values tools that lower barriers to chemical data access.
  Angle: MolBuilder's pure-Python approach aligns with his philosophy of accessible
  chemistry tools.

- **Cernak:** Reference the Cernak Lab's work on reaction discovery, high-throughput
  experimentation, or a specific paper on synthesis planning. He'll care about the
  reaction template coverage and scoring heuristics. Angle: MolBuilder could serve
  as a teaching companion for reaction informatics concepts.

- **Belford:** Reference OLCC or his work making chemistry education accessible
  online. He values open-source and community-driven tools. Angle: MolBuilder's
  Apache 2.0 license and Jupyter tutorials fit the open education model.

---

## Template B -- Python-in-Chemistry Educators

**Best for:** Faculty who already teach Python for chemistry/drug discovery.
Students can pip install immediately. Strong fit for computational labs.

**Targets:**

| Name | University / Org | Why they fit |
|------|-----------------|-------------|
| Charles Weiss | Wabash College | Teaches Python for chemists, authored "Scientific Computing for Chemists" |
| Pat Walters | OpenADMET / UMass | Drug discovery veteran, open-source advocate, runs Practical Cheminformatics blog |
| David Koes | U. Pittsburgh | Created 3Dmol.js, runs Koes Lab -- computational drug discovery and ML |

---

**Subject:** Pure-Python molecular engineering toolkit -- retrosynthesis to manufacturing

**Body:**

Hi [FIRST NAME],

I've been following [PERSONALIZE: "your Scientific Computing for Chemists textbook"
/ "your Practical Cheminformatics blog and OpenADMET work" / "3Dmol.js and the
Koes Lab's computational drug discovery tools"] -- it's been a real reference
for me.

I wanted to share MolBuilder, a Python toolkit I built that might complement what
you're teaching. The pitch: one `pip install` gives your students molecular
construction, SMILES parsing, 3D coordinate generation, retrosynthesis planning,
and process engineering -- with zero C++ dependencies.

Key details your students would care about:
- **Pure Python** (numpy/scipy/matplotlib) -- runs anywhere Python runs
- **270 purchasable starting materials** in the retrosynthesis knowledge base
- **24 functional group detectors**, 185 reaction templates
- **Jupyter tutorials** included: SMILES-to-3D, retrosynthesis, process engineering
- 1,280+ tests, Apache 2.0, Python 3.11/3.12/3.13

Install and try in 30 seconds:
```
pip install molbuilder
python -c "from molbuilder.smiles import parse; print(parse('c1ccccc1').name)"
```

Tutorials: https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/tutorials

Would love your thoughts -- especially on what would make it more useful for
teaching. Happy to add examples or adapt tutorials if there's interest.

Best,
Taylor Powell

GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
PyPI: https://pypi.org/project/molbuilder/

---

### Personalization Notes for Template B

- **Weiss:** Reference "Scientific Computing for Chemists" directly. He cares about
  pedagogy and progressive complexity. Angle: MolBuilder's tutorials follow a similar
  philosophy -- start with SMILES, build up to retrosynthesis, end at manufacturing.

- **Walters:** Reference his blog posts or OpenADMET. He's skeptical of hype and
  values honest engineering. Be upfront about limitations (rule-based, not ML).
  Angle: MolBuilder fills the process chemistry gap that drug discovery pipelines
  usually hand-wave past.

- **Koes:** Reference 3Dmol.js specifically -- MolBuilder's frontend uses it for
  visualization. He'll appreciate the technical details of the builtin 3D coordinate
  generator (distance geometry + force field). Angle: MolBuilder is a natural
  upstream data source for 3Dmol.js visualization workflows.

---

## Template C -- International Teaching Platforms

**Best for:** Groups that maintain notebook-based chemistry curricula. MolBuilder
can extend their existing platform with process chemistry coverage they don't have.

**Targets:**

| Name | University / Platform | Why they fit |
|------|----------------------|-------------|
| Andrea Volkamer | Charite Berlin | Created TeachOpenCADD -- open computational drug discovery teaching platform |
| Benjamin Morgan & Fiona Dickinson | U. Bath | Run pythoninchemistry.org -- Python teaching resources for chemistry |
| Erik Menke | UC Merced | Develops open computational chemistry course materials |

---

**Subject:** Process chemistry module for [PLATFORM NAME / your teaching notebooks]

**Body:**

Hi [PERSONALIZE: "Professor Volkamer" / "Ben and Fiona" / "Professor Menke"],

[PERSONALIZE: "TeachOpenCADD has been an incredible resource for computational drug
discovery education" / "pythoninchemistry.org has been a great reference for Python
chemistry teaching" / "Your open computational chemistry materials at UC Merced
have been helpful to me as both a learner and a developer"].

I built MolBuilder to cover an area I noticed most open teaching platforms don't
reach yet: **process chemistry** -- the bridge from a synthesis plan to actual
manufacturing decisions.

MolBuilder is a pure-Python toolkit (Apache 2.0, `pip install molbuilder`) that
covers:
- SMILES parsing and 3D coordinate generation
- Retrosynthesis with 185 reaction templates and beam search
- **Reactor selection, safety assessment, cost estimation, and scale-up** --
  the process engineering piece most curricula skip

It ships with Jupyter tutorials designed for classroom use:
https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/tutorials

I think it could complement [PERSONALIZE: "TeachOpenCADD's drug discovery
pipeline" / "your existing Python chemistry notebooks" / "your computational
chemistry course"] by adding process engineering as a downstream module -- showing
students what happens *after* you pick a synthetic route.

Would you be interested in taking a look? I'd be happy to develop additional
notebook content tailored to [PERSONALIZE: platform/course name].

Best,
Taylor Powell

GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
PyPI: https://pypi.org/project/molbuilder/

---

### Personalization Notes for Template C

- **Volkamer:** Reference TeachOpenCADD by name and a specific talktorial topic
  (e.g., T001 target identification or the ADMET modules). She values reproducibility
  and Jupyter-first workflows. Angle: MolBuilder adds the "what comes after lead
  optimization" step that TeachOpenCADD doesn't currently cover.

- **Morgan & Dickinson:** Reference pythoninchemistry.org directly and their
  approach of making computational methods accessible to experimentalists. They
  value simplicity and clear documentation. Angle: pure Python, no compilation,
  notebook-ready -- fits their philosophy perfectly.

- **Menke:** Reference his open course materials or computational chemistry teaching
  at UC Merced. He works at a smaller institution where open tools matter more.
  Angle: MolBuilder is free, well-tested, and ready to drop into an existing
  syllabus.

---

## Sending Checklist

Before sending each email:

1. [ ] Replace all `[PERSONALIZE]` and `[BRACKET]` placeholders
2. [ ] Verify the professor is still at the listed institution (check their faculty page)
3. [ ] Reference at least one specific paper, course, or project of theirs
4. [ ] Confirm tutorial links work: click through GitHub and PyPI links
5. [ ] Send from a professional email address
6. [ ] Track sends in a simple spreadsheet: name, date sent, response status
7. [ ] Follow up once after 7-10 days if no response (short, no pressure)

## Follow-Up Template (7-10 days, no response)

**Subject:** Re: [original subject]

**Body:**

Hi [NAME],

Just floating this back to the top of your inbox -- I know the semester is busy.
If MolBuilder isn't the right fit, no worries at all. If you do have 5 minutes to
look, the quickest way in is the first Jupyter tutorial:
https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/tutorials

Either way, thanks for the work you do in [PERSONALIZE: open chemistry education /
cheminformatics / computational drug discovery].

Best,
Taylor
