# IWOTA2023 slides

The [slides](iwota2023slides.pdf) for IWOTA2023 and the Julia codes that I used to produce pictures are available here.

This repository contains
* `README.md` : this file
* `.devcontainer.json` allows you to use this repository on GitHub Codespace and run the Julia code there
* `iwota2023slides.pdf` : slides I used in IWOTA2023
* `dysonBM.jl` simulates Dyson's Brownian motions of beta = 2
* `gff2D.jl` simulates Gaussian free field on a 2d square
* `twoDIsing_Dobrushin.jl` simulate the critical Ising model on a 2d square with Dobrushin boundary conditions (single interface)
* `twoDIsing_TwoInterfaces.jl` simulate the critical Ising model on a 2d square with boundary bonditions changing at the four corners (two interfaces)

I did not use `Plots` for Julia, but exported data as `.csv`.