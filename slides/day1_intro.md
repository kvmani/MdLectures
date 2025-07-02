
# Section 1: 1.1 – Welcome & Workshop Overview

## Slide: Title

- "Day 1: Intro to MD & LAMMPS for Metallurgy"
- Workshop Instructor, affiliation, date


## Slide: Objectives & Agenda

- Goals: grasp MD fundamentals, explore LAMMPS structure, run first demo
- Agenda:
  - short theory primer
  - demo overview
  - live LJ run
  - Ovito intro
  - Q&A

# Section 2: 1.2 – Molecular Dynamics Fundamentals

## Slide: What is MD?

- Definition: classical simulation of atoms using Newtonian mechanics
- Relevance to metallurgy and metals



## Slide: Equations of Motion

$$
\frac{d\mathbf{r}_i}{dt} = \mathbf{v}_i,\quad m_i\frac{d\mathbf{v}_i}{dt} = \mathbf{F}_i = -\nabla_i V(\{\mathbf{r}_j\})
$$

## Slide: Interatomic Potentials

- Lennard-Jones vs. EAM usage in metals
- Emphasize EAM for metallic bonding



# Section 3: 1.3 – LAMMPS Basics

## Slide: What is LAMMPS?

- Parallel MD code for materials simulation :contentReference[oaicite:1]{index=1}

## Slide: LAMMPS Workflow

- prepare input → run → analyze → visualize



## Slide: Anatomy of Input Script

```bash
# in.lj – Lennard-Jones demo
units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 10 0 10 0 10
create_box 1 box
create_atoms 1 box
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
fix 1 all nve
run 10000
```

# Section 4: 1.4 – Live Demo: LJ Melt in LAMMPS

## Slide: Setup & Goals

- LJ melt demo: simulate FCC melt, compute thermodynamic properties

## Slide: Input Script Overview

```yaml
include-code-files:
  - ../scripts/in.lj
```

## Slide: Run Commands

```bash
lmp -in in.lj
```

## Slide: Expected Output

- Thermo log sample




# Section 5: 1.5 – Visualization with Ovito

## Slide: Why Ovito?

- Supports LAMMPS dumps :contentReference[oaicite:2]{index=2}

## Slide: Loading Trajectories

will put ovito images later

## Slide: Basic Analysis in Ovito

- Displacement plots
- Coloring atoms by type or energy

# Section 6: 1.6 – Q&A & Day 2 Preview

## Slide: Recap

- MD basics, LJ demo, Ovito intro

## Slide: Day 2 Teaser

- crystal defects
- dislocations
- diffusion

