---
title: 'PyEAFE: A Python package to solve convection-diffusion-reaction equations using the edge-averaged finite element (EAFE) approximation'
tags:
  - Python
  - fenics
  - partial differential equations
  - convection-diffusion-reaction equations
authors:
  - name: Arthur Bousquet
    orcid: 0000-0003-2163-6473
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Maximilian Metti
    affiliation: "2"
  - name: Shuonan Wu
    affiliation: 3
affiliations:
 - name: Lake Forest College
   index: 1
 - name: Independent Researcher
   index: 2
 - name: Peking University 
   index: 3
date: 21 Frebruary 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx
aas-journal: Journal of Open Source Software
---

# Summary

This package is an implementation of the edge-averaged finite element (EAFE) approximation to linear convection-diffusion-reaction equations, which is use to stabilize discretizations for convection-domintated problems. To learn more about EAFE, refer to the EAFE section below.

Solvers for simulating charge-transport systems with an arbitrary number of charge-carrying species.
These systems are modeled by the Poisson-Nernst-Planck (PNP) equations with the possibility of coupling to the Navier-Stokes (NS) equation to simulate electrokinetic phenomena.
Our codes are based on the FEniCS package for discretizing the differential equations. This permits the use of any compatible linear algebra solver for solving the resulting linear system of equations.


# Statement of need

See [@Xu99amonotone] for more and see for [@BHMX18] applications. 

The equations are


$$
 -\nabla \cdot (\alpha(x) \nabla(u) + \beta(x) u) + R(x) u = f(x)
$$



# References
