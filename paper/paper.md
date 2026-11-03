---
title: 'UGME: A python package for non-Markovian protein folding dynamics'

tags:
  - Python
  - molecular dynamics
  - generalized master equations
  - non-Markovian dynamics
  - molecular kinetics

authors:
  - name: Anthony J. Dominic III
    orcid: 0000-0002-2240-0107
    affiliation: 1
    corresponding: true
  - name: Bertram Cham
    orcid: 0009-0005-7060-5157
    affiliation: 1
  - name: Andrés Montoya-Castillo
    orcid: 0000-0003-3037-3695
    affiliation: 1
    corresponding: true
affiliations:
  - name: Department of Chemistry, University of Colorado Boulder, Boulder, CO 80309, USA
    index: 1
date: 2 November 2026
bibliography: paper.bib
---

# Summary
Biomolecular processes, including protein folding, protein-protein association, and protein-DNA recognition, often occur on timescales beyond what individual atomistic molecular dynamics (MD) simulations can reach. To access these longer-timescale dynamics, practitioners commonly combine many parallelized, short trajectories into reduced kinetic models, such as Markov state models (MSMs) [@Bowman:2013; @Prinz:2011; @Husic:2018], that describe how the biomolecules interconvert between a handful of biomolecular configurations (e.g., folded and unfolded states of a protein). Central to MSM construction is the assumption that, at sufficiently long time lags, these projected dynamics become memoryless, such that the future evolution of the system depends only on its current state. Yet, projecting high-dimensional biomolecular dynamics onto a small set of states necessarily discards information about unresolved degrees of freedom, including solvent reorganization, side-chain rearrangements, and intermediate configurations. Because these unresolved motions can continue to influence subsequent dynamics, the resulting dynamics retain memory that is absent from a Markovian description. Here, we present `ugme`, a modular, open-source Python package for constructing and validating generalized master equations (GME) models that explicitly account for these non-Markovian (memory) effects.

![Schematic of the $\mathbf{U}$-GME approach for predicting long-time dynamics from short-time transition probability matrices (TPMs). **Top:** A time series of 3$\times$3 TPMs, $\mathbf{T}(t)$, represented as heatmaps, beginning with the identity matrix $\mathbf{T}(0) = \mathbf{I}$ and extending to the final available time $N$ with timestep $\Delta t$. By applying the `ugme` software, one can predict the TPM dynamics for times beyond $N$. The black-outlined matrix elements trace the conditional probability of remaining in state 3 shown in the bottom-left panel. **Bottom left:** Time evolution of $\mathbf{T}_{33}(t)$,
 illustrating that the dynamics are still relaxing over the interval defined by the averaging onset $t_r$ and cutoff $\tau_R$. **Bottom right:** Corresponding evolution of the propagator element $\mathbf{U}_{33}(t)$ (in light purple), which approaches a stationary value despite the continued evolution of the TPM. The propagator is averaged over the shaded interval $[t_r, \tau_R]$ and subsequently held constant for $t > \tau_R$, enabling prediction of dynamics beyond the available TPM dataset.](schematic.png)
 
# Statement of need
`ugme` is designed to extend existing molecular kinetics workflows to systems whose projected dynamics retain memory. Current molecular simulation analysis software packages provide well-established tools for trajectory featurization, dimensionality reduction, clustering, and transition probability matrix (TPM) estimation. With a time series of TPMs in hand, existing software packages provide robust methods for estimating, validating, and analyzing MSMs. However, once the dynamics have been projected onto a reduced set of states, accounting explicitly for emergent memory effects generally requires specialized implementations not integrated into standard analysis pipelines. This motivates the development of a tool that can operate directly on TPMs produced by existing workflows.

`ugme` addresses this need by implementing the noise-resilient time-local GME framework, the $\mathbf{U}$-GME introduced in @Dominic:2023a. As we show in Fig. 1, given a time series of TPMs from simulation workflows, analytical models, or experimental sources, the package computes the time-dependent propagator that captures the dynamical memory effects. `ugme` then identifies when these propagators have converged and uses their long-time behavior to predict dynamics beyond the total time spanned by the input dataset. The package also provides diagnostics for assessing convergence and testing the sensitivity of predicted kinetics to choices such as the memory relaxation time and propagator averaging window.

Because `ugme` operates directly on TPMs, it complements existing MD analysis packages such as MSMBuilder, PyEMMA, and deeptime. Users can continue to employ established tools for featurization, clustering, and TPM estimation, while applying `ugme` alongside MSMs to improve temporal resolution and quantify memory effects from unresolved degrees of freedom, or instead of MSMs when long lag times or underconverged statistics preclude their use.

# State of the field
Established packages for molecular kinetic modeling primarily focus on Markovian descriptions. MSMBuilder [@Harrigan:2017], PyEMMA [@Scherer:2015], and deeptime [@Hoffmann:2022] provide workflows for estimating Markov models from simulation data, while msmhelper [@Nagel:2023] provides lightweight tools for estimating, validating, and analyzing MSMs from discrete trajectories. These packages are well-suited to systems for which a sufficiently Markovian state representation and lag time can be identified. However, in @Dominic:2023a, we demonstrate that a Markovian description is insufficient for capturing mRNA recognition dynamics by the human Argonaute complex.

Recent developments have extended molecular kinetic modeling beyond the Markovian approximation using GMEs which retain the influence of unresolved degrees of freedom on the projected dynamics [@Zwanzig:1960; @Mori:1965] Memory kernel formulations encode this influence through a time-convolution integral, whereas time-local formulations describe it through an explicitly time-dependent propagator or generator, akin to a time-dependent rate matrix. While the two representations are formally similar, their numerical implementations differ substantially. For example, the integrative GME (IGME) uses time-integrated memory kernels to construct numerically stable models of long-time biomolecular dynamics [@Cao:2023], for which a tutorial and accompanying code exists [@Wu:2024] These developments demonstrate the growing utility of GMEs for biomolecular dynamics and motivate complementary software implementations based on alternative formulations of the GME.

`ugme` implements a time-local construction that requires only a sequence of TPMs. This formulation avoids direct evaluation of projected forces, numerical differentiation of correlation functions, or explicit convolution of a memory kernel. The scope of `ugme` therefore differs from existing MSM packages. While existing tools provide machinery for constructing and analyzing Markovian kinetic models, `ugme` begins from the resulting TPMs and analyzes the memory retained in their dynamics. Keeping these steps separate allows `ugme` to remain compatible with TPMs produced by existing software while providing a workflow centered specifically on the estimation, validation, and propagation of non-Markovian kinetic models.

# Software design

The primary input to `ugme` is a time-indexed sequence of square TPMs, $\{ \mathbf{T}(t_0),\mathbf{T}(t_1),\ldots,\mathbf{T}(t_N) \}$, where each $\mathbf{T}(t_n)$ contains the conditional transition probabilities between the states of the reduced model at time $t_n = n \Delta t$, where $\Delta t$ is the time interval between successive TPMs (i.e., the temporal resolution of the reference dataset) and $t_N$ is the final timestep of the input TPM dataset. Crucially, their time evolution evolves according to $\mathbf{T}(t_{n+1}) = \mathbf{U}(t_n) \mathbf{T}(t_n)$, where the time-dependent propagator $\mathbf{U}(t_n)$ provides a direct measure of the relaxation of memory in the reduced dynamics. The memory in $\mathbf{U}(t_n)$ decays and non-Markovian effects become unimportant when its matrix elements reach a steady state. As such, the memory relaxation timescale $\tau_R \equiv t_R = R \Delta t$ is the time at which $\mathbf{U}(t_n)$ plateaus. For times beyond $\tau_R$, we replace the propagator with its stabilized value $\mathbf{U}_\infty (\Delta t)  \equiv \mathbf{U}(\tau_R)$. Beyond $\tau_R$, the dynamics evolve according to
$$
\mathbf{T}(\tau_R + n\Delta t) = \mathbf{U}_\infty^n(\Delta t)\mathbf{T}(\tau_R).
$$

The package is organized around the `UGME` model object, which stores the TPM time series, temporal resolution, and quantities generated during model construction. From the TPMs, `ugme` computes the corresponding time-local propagator $\mathbf{U}$ by directly solving the time-local form for $\mathbf{U}(t_n) = \mathbf{T}(t_{n+1}) [\mathbf{T}(t_{n})]^{-1}$ and tracks its time evolution. The plateauing of $\mathbf{U}(t)$ quantifies the memory's relaxation time in the reduced dynamics. For sufficiently converged data, `ugme` identifies the memory-relaxation timescale $\tau_R$ using the RMSE-based convergence criterion we established in @Dominic:2023a.

TPMs estimated from finitely-sampled trajectory data can introduce statistical error into the late-time propagator. To reduce the propagator's sensitivity to this noise, `ugme` supports time-averaging of the propagator $\mathbf{U}$. We introduce tools to compute the onset ($t_r$) and offset ($\tau_R$) averaging parameters and the averaged propagator that we define as
$$ \langle \mathbf{U}_\infty(\Delta t) \rangle  = \frac{1}{R - r}\sum_{n = r}^{R-1} \mathbf{U}(t_n).$$
To suppress noise and propagate the dynamics, we insert the averaged propagator into our GME to obtain $$ \mathbf{T}(\tau_R + n\Delta t) = \langle \mathbf{U}_\infty(\Delta t) \rangle ^n \mathbf{T}(\tau_R).$$

In a typical workflow, a user supplies a TPM time series, constructs the corresponding `UGME` model, determines $\tau_R$, selects an averaging window when needed, and propagates to access the long-time dynamics. This organization keeps the model diagnostics accessible while maintaining a compact interface for routine calculations. These stages are exposed separately so that users can inspect the quantities entering the model and assess if and how numerical choices affect the resulting kinetics.

The repository includes source code, API documentation, example NumPy datasets, and a notebook demonstrating TPM loading, propagator construction, reconstruction, relaxation-time scans, averaging-window scans, and long-time prediction. These materials provide a reproducible starting point for applying the method to new biomolecular systems and for extending the implementation to additional GME estimators.

# Research impact statement

We developed the methods implemented in `ugme` to extract long-time biomolecular kinetics from short-time simulation data. In @Dominic:2023a, we applied the $\mathbf{U}$-GME to alanine dipeptide, the FiP35 WW Domain, and the human Argonaute complex, showing that the method can recover long-time relaxation dynamics from substantially shorter input trajectories than required by the conventional MSM treatment. We subsequently discussed the broader role of memory-enriched models for biomolecular simulations in @Dominic:2023b.

Since then, researchers have adopted their own time-local formulation. For example, @Liu:2026 has adapted the time-local GME formulation to estimate discrete time-local GMEs from experimental and simulation data by constructing physically constrained maximum likelihood estimators and applied these methods to experimental datasets corresponding to single-molecule Förster resonance energy transfer (FRET) in an ion channel, experimental nanoparticle tracking of a processive molecular motor, and the simulated protein folding dynamics of HP35. Beyond biophysics, the time-local GME has extended to polymer dynamics,[@Bement:2025; @Sinclair:2025] transport,[@Bhattacharyya:2025a; @Bhattacharyya:2025b] and spectroscopy [@Sayer:2026] Collectively, these studies illustrate that short-lived memory can be exploited across a range of dynamical problems where projecting onto a reduced set of variables introduces memory effects.

By separating GME analysis from TPM estimation, `ugme` makes these explicit approaches to capturing memory effects readily available for analyzing molecular dynamics. We anticipate that this interface will facilitate applications of time-local GME models to new biomolecular systems while also providing a foundation for testing and extending memory-enriched kinetic modeling methods.

# AI usage disclosure

A significant portion of this code was developed from 2021-2023 without the use of AI tools. AI tools were used to proofread code, proofread the publication, and build initial docstrings and comments for portions of the code that previous implementations of this code did not have.

# Data Availability
Our materials are freely available at the following GitHub repository: https://github.com/ajdominic/ugme.github.io/tree/main.

# Acknowledgments

A.J.D., B.C., and A.M.C. acknowledge support by an Early Career Award in the CPIMS program in the Chemical Sciences, Geosciences, and Biosciences Division of the Office of Basic Energy Sciences of the U.S. Department of Energy under Award DE-SC0024154. A.J.D. also acknowledges partial support from the NIH/CU Molecular Biophysics Program and the NIH Biophysics Training Grant T32 GM145437. This work was supported by an Alfred P. Sloan Research Fellowship and the Camille and Henry Dreyfus Foundation through a Camille Dreyfus Teacher-Scholar Award. This work utilized the Alpine high performance computing resource at the University of Colorado Boulder. Alpine is jointly funded by the University of Colorado Boulder, the University of Colorado Anschutz, Colorado State University, and the National Science Foundation (Award No. 2201538).

# References
