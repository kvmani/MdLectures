# Hands-on Molecular Dynamics with LAMMPS — Slide Deck & Workshop Repo

This repository hosts the **Quarto-based slide decks, example input files, and GitHub Actions workflows** for a two-day online workshop on classical molecular dynamics (MD) using **LAMMPS** (Large-scale Atomic/Molecular Massively Parallel Simulator).  
The project is designed so that:

* you can **clone / fork** and immediately render the slides to HTML (Reveal.js), PDF (Beamer), or PPTX with a single `quarto render`,  
* GitHub Actions automatically rebuilds and publishes the latest slides to **GitHub Pages** (`docs/`), so the live deck is always available at  
  `https://kvmani.github.io/MdLectures/slides/intro.html`,  
* attendees can download or run the accompanying LAMMPS input scripts directly from the repo.

---

## 📂 Directory layout

```text
MdLectures/
├── _quarto.yml          # project config → outputs to docs/
├── images/              # all slide figures & QR codes
│   ├── lammps_logo_placeholder.png
│   ├── md_integration_loop_placeholder.png
│   └── ...
├── slides/              # each lecture = one .qmd file
│   └── intro.qmd
├── docs/                # rendered site (auto-generated by CI)
├── .github/
│   └── workflows/
│       ├── build.yml    # builds HTML/PDF/PPTX on every push
│       └── publish.yml  # (optional) deploys docs/ to GitHub Pages
└── .gitignore

```

## 🚀 Quick start

### Clone and install Quarto
```bash
git clone https://github.com/kvmani/MdLectures.git
cd MdLectures
quarto check
```

### Render slides locally
```bash
quarto render slides/intro.qmd   # HTML → docs/slides/intro.html`

```

### Preview website
```bash
quarto preview
```

## 🔧 Tools you'll need

- **LAMMPS** – compiled with MPI support for running the example scripts.
- **OVITO** – for visualising dump files generated during demos.
- **Python 3** – used for small utilities and compatible with Quarto.
- **Git** – clone this repository and track your own input scripts.
- **Quarto** – to render the slides to HTML, PDF or PPTX formats.

## 📚 How to follow the lecture series

1. Clone the repository and install the tools listed above.
2. Run `quarto render slides/day1_intro.qmd` to generate the first lecture in your preferred format.
3. Explore the `scripts/` folder and run the accompanying LAMMPS inputs, e.g. `lmp -in scripts/in.lj`.
4. Visualise the output dumps using Ovito to reinforce concepts covered in the slides.
5. Proceed through each `slides/*.qmd` file in numerical order—each corresponds to a workshop day.


## Acknowledgements

* LAMMPS development team at Sandia National Labs and the user community.

* Quarto open-source authors for the flexible publishing toolchain.

* GitHub for CI/CD & Pages hosting.
