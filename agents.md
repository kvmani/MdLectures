# Agents Overview

This file defines repository structure, purposes, and requirements. It guides Codex as the "agent" responsible for generating lecture materials, slides, tutorial code, and build pipelines.

---

## 🧠 Purpose

- **Primary goal**: Generate materials (slides, scripts, tutorials) for 2×2‑hour hands‑on LAMMPS-based Molecular Dynamics workshops focused on metallurgical problems (dislocations, defects, diffusion).
- **Target audience**: PhD scholars and junior faculty in materials science.
- **Agent tasks**:
  - Create slide content in Quarto (.qmd) covering theory, examples, and live demos.
  - Populate `scripts/` folder with LAMMPS scripts, richly commented for clarity.
  - Embed code snippets in slides using standard code font (```).
  - Manage images and assets referenced in slides.
  - Where ever images are to be inserted, use place holder .svg files which will be later changed to actual images.
  - avoid committing any binary files (like .png, .pdf, .pptx etc at all costs.)

---

## 📁 Repository Structure
MdLectures/
├── agents.md # This file
├── slides/ # All Quarto slides (.qmd)
├── images/ # All images used in slides
├── scripts/ # LAMMPS input scripts for demos
├── _quarto.yml # Quarto project config
├── requirements.txt # Build dependencies
└── README.md # Repo overview & workshop logistics


---

## Folder Responsibilities

### `slides/`
- Contains `.qmd` files for lectures.
- Use clear YAML front matter: specify title, author, formats (revealjs, pptx, beamer, html).
- Embed theoretical explanations, code chunks (```), and slide directives.
- Refer to images located in `images/`.

### `images/`
- Stores all assets referred from slides.
- Use descriptive filenames (`dislocation_twin.png`, `diffusion_plot.svg`, etc.).

### `scripts/`
- Contains LAMMPS input scripts (`.in`, `.lmp`) used in demos.
- **Mandatory**: each script includes comments explaining:
  - Purpose of simulation
  - Key commands and parameters
  - Expected output / how to visualize

---

## ⚙️ Build Pipeline

### Quarto config (`_quarto.yml`)
- Declare project formats:
  - `revealjs` (HTML slides)
  - `pptx` (PowerPoint)
  - `beamer` (LaTeX/PDF)
  - `html` (web)
- Ensure `self-contained: true` for revealjs for easy sharing :contentReference[oaicite:1]{index=1}.
- Specify monospace `monofont` for code snippets in all formats :contentReference[oaicite:2]{index=2}.

### `requirements.txt`
- List build tools and dependencies:
  - `quarto`
  - `pandoc`
  - LaTeX engine: `tinytex`
  - (Optional) `revealjs`, `node.js` for advanced HTML
  
---

## ✅ Best Practices & Agent Instructions

1. **Modular content**: Separate slide chapters by topic (e.g. `day1_intro.qmd`, `day2_diffusion.qmd`).
2. **Include code files**:
   ```yaml
   include-code-files:
     - ../scripts/dislocation_demo.in

Each .in script begins with a header comment block:

# dislocation_demo.in
# Purpose: simulate edge dislocation in FCC lattice
# Steps:
#   - create geometry
#   - apply deformation
#   - compute stress/energy

Dependencies check: Agent should verify required executables are installed before generating materials.

Output formats: On render, agent ensures .pptx, .pdf, .html are generated and placed in _site/ or docs/.

Prioritize content quality over styling—clear MD slides with well‑commented code.

All maths outputs must be verified to ensure they got compiled well (look for tex logs and errors and correct them when necessary)

Enforce folder structure: slides/, scripts/, images/.

Use Quarto config to produce multiple formats.

Embed code in slides, auto-include script files, and use code font.

Scripts in scripts/ must be self-explanatory through comments.

Ensure build reproducibility with requirements.txt.