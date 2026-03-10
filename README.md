# 3-Component Radioactive Decay Chain Simulator

Simulates a 3-isotope decay chain (A → B → C) using both an analytical solution and a forward Euler numerical approximation. Generates population plots and a convergence analysis.

## Requirements

- [Julia](https://julialang.org/downloads/) (tested on 1.x)
- Julia packages: `Plots`

Install the required package by running in the Julia REPL:

```julia
using Pkg
Pkg.add("Plots")
```

## Running with the Default Input

The repo includes a sample `input.txt`. To run:

```bash
julia main.jl
```

This will produce:
- `output.txt` — time-series table of NA, NB, NC populations
- `plot1_convergence.png` — NB numerical convergence at multiple step sizes
- `plot2_populations.png` — all populations (numerical, Δt = 0.125 h)
- `plot2b_populations_analytical.png` — all populations (analytical)
- `plot3_tmax.png` — convergence of the time at which NB peaks

## Using Your Own Input

Edit `input.txt` in the project directory. The format is one parameter per line:

```
hl_A    <half-life of isotope A in hours>
hl_B    <half-life of isotope B in hours>
Na0     <initial population of A>
Nb0     <initial population of B>
Nc0     <initial population of C>
tfinal  <simulation end time in hours>
```

Example (the default):

```
hl_A 1.97
hl_B 9.91
Na0 100.0
Nb0 0.0
Nc0 0.0
tfinal 50.0
```

Then run `julia main.jl` as above.

## File Overview

| File | Description |
|---|---|
| `main.jl` | Entry point: reads input, runs simulations, writes output and plots |
| `projFunctions.jl` | Analytical solutions, numerical integrator, and helper functions |
| `input.txt` | Simulation parameters |
| `output.txt` | Generated time-series output |
