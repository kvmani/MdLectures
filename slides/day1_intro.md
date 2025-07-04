---
title: "Molecular Dynamics with LAMMPS – Day 1"
author: "KV Mani Krishna"
date: 2025-07-02
---

## What is Molecular Dynamics? 🧬
- Simulates the time evolution of **N** interacting particles by integrating Newton’s equations (F = ma).  
- Forces come from an interatomic **potential/force-field** that sets the simulation’s fidelity.  
- Reaches picoseconds → microseconds and nanometres → microns on modern hardware.
![MD then vs now](../images/md_history_vs_today.png)

**Take-away** MD is a digital time-lapse microscope for atoms.

---

## The MD Algorithm at a Glance ⚙️
- **Velocity-Verlet loop**: ½-kick → drift → (re)build neighbours → forces → ½-kick.
- **Neighbour lists + domain decomposition** keep force evaluation O(N).
- Typical timestep ≈ 1 fs; millions of steps deliver nanosecond trajectories.

**Pro-tip** Let physics—not code—set the timestep.

---

## Problems MD Can Tackle 🛠️
- Point & line defects: vacancies, dislocations, radiation damage.
- Phase transitions: melting, solidification, diffusion-driven growth.
- Mechanical & thermal properties: elastic moduli, crack propagation, conductivity.

**Key takeaway** If atoms move or rearrange ≤ µm & ≤ µs, MD can probably watch it.

---

## LAMMPS Philosophy & Structure 🏗️
- Open-source C++ with modular **styles** (pair, fix, compute, …) for plug-in physics.
- Continuous-release; extra packages can load at runtime as shared libs.
- Scales from laptops to exascale GPUs via domain decomposition + MPI.
- Usable standalone or as a callable library/Python module.

![Lammps Logo](../images/lammps_logo_placeholder.png)


**Take-away** LAMMPS is Lego® for materials simulation.

---

## LAMMPS Scripting 101 ✍️
- Plain-text scripts; supports variables, math, loops, conditionals, shell cmds.
- Flow: `units` → build atoms → potentials → `fix`es → `run` → output.
- 1000+ commands; docs & ~35 example folders live at lammps.org.

**Pro-tip** Treat every *.in* file like code—version-control and iterate.

---

## Anatomy of a LAMMPS Input 🩺
- **Setup**: `units`, `dimension`, `boundary`, `atom_style`, `neighbor`.
- **Geometry**: `lattice`, `region`, `create_box|read_data`, `create_atoms`.
- **Interactions & groups**: `pair_style`, `pair_coeff`, `group`, `region`.
- **Dynamics**: `mass`, `velocity`, `fix nve/nvt/npt`, thermostats/constraints.
- **Output & run**: `timestep`, `thermo`, `dump`, `run`, `write_restart`.

**Take-away** Five sections, infinite possibilities.

---

---

## Atom Styles & Potentials in LAMMPS  
.key-takeaway[
LAMMPS provides flexible `atom_style` options to match the physics of your system, and a wide range of interatomic potentials suited for metals, semiconductors, and molecules.
]

- The `atom_style` determines what per-atom properties are stored (e.g., charge, bonds, spin).
- The `pair_style` defines how atoms interact — ranging from simple to complex potentials.
- Matching the correct pair of `atom_style` and `pair_style` is crucial for meaningful simulations.

---

## Common `atom_style` Options  


| `atom_style`   | Description                                 |
|----------------|---------------------------------------------|
| `atomic`       | For simple atoms (metals, no charges/bonds) |
| `charge`       | Adds charge for ionic systems               |
| `molecular`    | Includes bonds, angles, dihedrals           |
| `full`         | `molecular` + charge (e.g., biomolecules)   |
| `sphere`       | Adds size and rotation (colloids)           |
| `ellipsoid`    | Non-spherical particles                     |
| `hybrid`       | Combine multiple styles                     |

.key-takeaway[
Use minimal styles like `atomic` for metals, and extended ones like `full` or `molecular` for complex molecules or charged species.
]
---

## Common Potentials and Their Use Cases  


| Potential       | `pair_style`     | Suitable For                  | Notes                                    |
|------------------|------------------|-------------------------------|------------------------------------------|
| Lennard-Jones    | `lj/cut`         | Noble gases, coarse-grained   | Simple 12-6 interaction                   |
| EAM              | `eam`, `eam/alloy` | FCC/BCC metals (e.g., Cu, Fe) | Many-body; ideal for pure metals         |
| Tersoff          | `tersoff`        | Covalent solids (Si, C)       | Bond-order; angular dependence           |
| Stillinger-Weber | `sw`             | Semiconductors (Si, Ge)       | Three-body; angular term                 |
| MEAM             | `meam`           | Alloys, metals                | Angular EAM; supports complex alloys     |
| ReaxFF           | `reax/c`         | Reactive systems, combustion  | Break/form bonds; needs charge equil.    |
| Morse            | `morse`          | Diatomic bonded atoms         | Smooth potential for simple molecules    |
| Buckingham       | `buck`           | Ionic crystals, oxides        | Used with Coulombic terms                |

.key-takeaway[
Choose potentials that reflect bonding physics: EAM for metals, Tersoff for covalent systems, and ReaxFF for reactive systems.
]

## Example – LJ Melting 🔥
```lmp
# in.lj_melt — LJ solid → liquid
units           lj
dimension       3
boundary        p p p
atom_style      atomic

lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass 1 1.0


pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5

velocity        all create 0.3 87287
fix             1 all nve
thermo          100
thermo_style    custom step temp pe etotal press

dump            1 all atom 200 dump.lj
timestep        0.005
run             20000
```

## Advanced Example: Lattice a(T) – V1→V4 📈

### V1 – single-T NPT & print
```lmp
variable T equal 300
fix 1 all npt temp ${T} ${T} 0.1 iso 0 0 1.0
variable a equal lx/10
thermo_style custom step temp v_a
```
### V2 – temperature loop
```lmp
variable T loop 300 600 900 1200
label loopT
  fix 1 all npt temp ${T} ${T} 0.1 iso 0 0 1.0
  run 20000
  print "${T} $(v_a)" append latparam.dat
  unfix 1
next T
jump SELF loopT
```
### V3 – on-the-fly compute
```lmp
compute lat all reduce ave c_myCell[1]
variable a equal c_lat/10
fix 2 all ave/time 100 10 1000 v_a file lat_vs_t.dat
```
**Take-away**  Evolve scripts incrementally—variables & loops unlock automation.


# Section 2.1 – Problem Statement & V1 Script

## Fixed‑parameter script (`in.script_V1`)

```lammps
# in.script_V1 — lattice‑parameter at 300 K, 6×6×6 FCC Al cell
units           metal
atom_style      atomic
boundary        p p p
lattice         fcc 4.05
region          box block 0 6 0 6 0 6
create_box      1 box
create_atoms    1 box

pair_style      eam
pair_coeff      * * Al_u3.eam

timestep        0.002
velocity        all create 300.0 12345

fix             1 all npt temp 300 300 0.1 iso 0 0 1.0
thermo_style    custom step temp press vol lx
thermo          500
run             10000          # 20 ps
```

**Run command**

```bash
lmp -in in.script_V1
```

> **Pro‑tip:** hard‑wiring parameters means four manual edits just to change the temperature!

---

# Introducing Variables

## Script with variables (`in.script_V2`)

```lammps
# --- user‑defined variables ---
variable  a0    equal 4.05      # Å
variable  T     equal 300       # K
variable  nx    equal 6
variable  steps equal 10000     # 0.002 ps × 10 000 = 20 ps
# ------------------------------

units           metal
atom_style      atomic
boundary        p p p
lattice         fcc ${a0}
region          box block 0 ${nx} 0 ${nx} 0 ${nx}
create_box      1 box
create_atoms    1 box

pair_style      eam
pair_coeff      * * Al_u3.eam

timestep        0.002
velocity        all create ${T} 12345

fix             1 all npt temp ${T} ${T} 0.1 iso 0 0 1.0
thermo_style    custom step temp lx
thermo          500
run             ${steps}
```

**Run with overrides**

```bash
lmp -var T 500 -var nx 8 -in in.script_V2
```

> **Pro‑tip:** `-var name value` on the command line lets you sweep parameters without editing the file.

---

# Loop over Temperature

## Temperature sweep (`in.script_V3`)

```lammps
# temperature list
variable  Tlist index  200 300 400 500 600
variable  nx    equal 6
variable  a0    equal 4.05

label loopT
variable T equal ${Tlist}

lattice         fcc ${a0}
region          box block 0 ${nx} 0 ${nx} 0 ${nx}
create_box      1 box
create_atoms    1 box
pair_style      eam
pair_coeff      * * Al_u3.eam
timestep        0.002
velocity        all create ${T} 12345
fix             1 all npt temp ${T} ${T} 0.1 iso 0 0 1.0

log             log_T${T}.lammps
run             10000
undump          all
clear
next Tlist
jump SELF loopT
```

**Execution**

```bash
lmp -in in.script_V3
```

> **Pro‑tip:** use `clear` before `next` to avoid “system already exists” errors in loops.

---

# Self‑describing Filenames

## String variables (`in.script_V4`)

```lammps
variable  Tlist index 300 400 500
variable  nx    equal 6
variable  a0    equal 4.05

label loopT
variable T equal ${Tlist}

log    logs/log_T${T}.lammps
variable dumpfile string dumps/dump_T${T}_N${nx}.lammpstrj
dump   1 all custom 200 ${dumpfile} id type x y z

# ... simulation setup ...
run 10000
undump 1
clear
next Tlist
jump SELF loopT
```

File names now carry temperature and cell size (`dump_T300_N6.lammpstrj`).

---

# Complete Sweep & CSV Output

## Final script (`in.script_final`)

Key additions: loop over cell sizes and write lattice parameter to CSV.

```lammps
variable Tlist index 300 400 500
variable nlist index 4 6 8
variable a0 equal 4.05

label loopN
variable nx equal ${nlist}

label loopT
variable T equal ${Tlist}

# build cell
lattice fcc ${a0}
region box block 0 ${nx} 0 ${nx} 0 ${nx}
create_box 1 box
create_atoms 1 box
pair_style eam
pair_coeff * * Al_u3.eam
timestep 0.002
velocity all create ${T} 12345
fix 1 all npt temp ${T} ${T} 0.1 iso 0 0 1.0

# lattice parameter output
variable lat equal lx/${nx}
fix avg all ave/time 100 1 100 v_lat file results/a_vs_T.csv mode append

run 10000
unfix avg
undump all
clear

next Tlist
jump SELF loopT

label nextN
next nlist
jump SELF loopN
```

Run all sweeps:

```bash
lmp -in in.script_final
```

---

# Recap & Preview

- Variables and loops turn a rigid script into a reusable driver.  
- Self‑describing filenames eliminate confusion during post‑processing.  
- Next: defect formation energies, dislocation mobility, diffusion studies.

---
### Script V1 – Hard-wired
```lammps
# in.script_V1 -- Hard‑wired single‑temperature run
# Purpose: demonstrate a minimal LAMMPS input that measures lattice parameter
# at 300 K for a 6×6×6 FCC Al cell.  No variables, no loops.

units           metal
atom_style      atomic
boundary        p p p

# --- system definition -------------------------------------------------------
lattice         fcc 4.05                      # a0 fixed
region          box block 0 6 0 6 0 6
create_box      1 box
create_atoms    1 box

# --- interactions ------------------------------------------------------------
pair_style      eam
pair_coeff      * * Al_jnp.eam                # <‑‑ changed to Al_jnp.eam

# --- run setup ---------------------------------------------------------------
timestep        0.002         # ps
velocity        all create 300 12345
fix             1 all npt temp 300 300 0.1 iso 0 0 1.0
thermo_style    custom step temp lx
thermo          500

# --- run ---------------------------------------------------------------------
run 10000       # 20 ps
```
[open original](../scripts/in.script_V1)

### Script V2 – Variables …
```lammps
# in.script_V2 -- Introduce user variables
# Changes from V1:
#  1. a0, T, nx, and steps promoted to variables so they can be overridden
#     via '-var name value' on the command line.
#  2. Region size and thermostat now reference these variables.
#  3. Still single run, single output; no loops yet.

# --- user‑defined variables (default values) ---------------------------------
variable  a0    equal 4.05      # Å
variable  T     string 300       # K  (string style so cmd‑line -var overrides)
variable  nx    string 6
variable  steps string 10000     # timesteps
# -----------------------------------------------------------------------------

units           metal
atom_style      atomic
boundary        p p p

lattice         fcc ${a0}
region          box block 0 ${nx} 0 ${nx} 0 ${nx}
create_box      1 box
create_atoms    1 box

pair_style      eam
pair_coeff      * * Al_jnp.eam

timestep        0.002
velocity        all create ${T} 12345

fix             1 all npt temp ${T} ${T} 0.1 iso 0 0 1.0
thermo_style    custom step temp lx
thermo          500
run             ${steps}

print "INFO: Completed simulation at T=${T} K for ${steps} steps"
```
[open original](../scripts/in.script_V2)

### Script V3 – Loop over Temperature
```lammps
# in.script_V3 -- Loop over temperature list, per‑T files
# New features versus V2:
#  * variable Tlist index ...; loop with next/jump
#  * self‑describing log/dump filenames containing T
#  * clear/undump so that loop restarts clean
#  * still single cell size

variable a0     equal 4.05
variable nx     string 6
variable steps  string 2000
variable Tlist  index 200 300 400 500 600

label loopT
variable T equal ${Tlist}

log             logs/log_T${T}.lammps
variable dumpfile string dumps/dump_T${T}_N${nx}.lammpstrj

units           metal
atom_style      atomic
boundary        p p p
lattice         fcc ${a0}
region          box block 0 ${nx} 0 ${nx} 0 ${nx}
create_box      1 box
create_atoms    1 box
pair_style      eam
pair_coeff      * * Al_jnp.eam

### for calcualting the lattice parameter
variable lat equal lx/${nx} 

timestep        0.002
velocity        all create ${T} 12345
fix             1 all npt temp ${T} ${T} 0.1 iso 0 0 1.0
dump            1 all custom 200 ${dumpfile} id type x y z
thermo_style    custom step temp lx
thermo          500
run             ${steps}

#print "INFO: Completed simulation at T=${T} K, lattice #parameter stored in ${dumpfile}" append yes

print "RESULT: T=${T} K nx=${nx} a_eq=${lat} dump=${dumpfile}" append results/results.txt  screen yes

undump 1
clear
next Tlist
jump SELF loopT
```
[open original](../scripts/in.script_V3)

### Script V4 – Size & Temperature Sweep
```lammps
# in.script_final — size & temperature sweep, NO if-statements

variable a0     equal 4.05
variable nx     index 4 6 8                  # outer loop list
variable steps  equal 2000

# -------- outer loop over nx -----------------------------------------
label loopNx

	# define temperature list fresh for each nx
	variable T delete
	variable T index 300 400 500 600             # inner loop list

	# -------- inner loop over T ------------------------------------------
	label loopT

		log       logs/log_T${T}_N${nx}.lammps
		variable  dumpfile string dumps/dump_T${T}_N${nx}.lammpstrj

		units     metal
		atom_style atomic
		boundary  p p p
		lattice   fcc ${a0}
		region    box block 0 ${nx} 0 ${nx} 0 ${nx}
		create_box 1 box
		create_atoms 1 box

		pair_style eam
		pair_coeff * * Al_jnp.eam

		# minimisation
		min_style cg
		minimize 1e-6 1e-8 1000 10000

		# NPT equilibration
		timestep   0.002
		velocity   all create ${T} 12345
		fix        1 all npt temp ${T} ${T} 0.1 iso 0 0 1.0

		dump       1 all custom 200 ${dumpfile} id type x y z
		variable   lat equal lx/${nx}
		fix        avg all ave/time 100 1 100 v_lat append results/a_vs_T.csv

		thermo_style custom step temp press v_lat
		thermo     500
		run        ${steps}

		print "RESULT: T=${T} K nx=${nx} a_eq=${lat} Ang dump=${dumpfile}" append results/results.txt screen yes

		undump 1
		unfix 1
		unfix avg
		clear

	next T
	jump SELF loopT            # go back unless temperature list is finished

# -------- advance cell size ------------------------------------------
next nx
jump SELF loopNx            # go back unless nx list is finished
```
[open original](../scripts/in.script_V4)
