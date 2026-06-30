# m-th King Method for Vortex PAL Modeling

This repository contains MATLAB codes for modeling vortex-based parametric array loudspeakers (PALs) using an m-th order King-integral framework. The codes include ultrasonic-field calculation, audible-field calculation, convergence analysis, direct-integration-method (DIM) comparison, local-effect correction, and computational-cost benchmarking.

## Overview

Conventional King-integral formulations are mainly suited to axisymmetric PAL fields. For vortex-based PALs, nonzero topological charges introduce azimuthal phase dependence and lead to non-axisymmetric fields. This repository provides MATLAB implementations for calculating such fields using an m-th order King/FHT-based framework and for comparing the results with direct integration where applicable.

The code is intended for research use and for reproducing numerical results related to vortex PAL modeling.

## Main Scripts

- `JASA_Calc_convergence_ultra.m`  
  Convergence analysis for ultrasonic fields.

- `JASA_Calc_convergence_audio.m`  
  Convergence analysis for audible fields.

- `JASA_KingDIM_lineCompare.m`  
  Line comparison between the proposed King/FHT method and the direct integration method, including local-effect correction.

- `JASA_KingDIM_ComplexityCompare.m`  
  Timing, memory, and sample-count benchmark for the proposed method and DIM.

- `JASA_AudibleFieldSweepm.m`  
  King-only audible-field sweep for different modal-order combinations with King local-effect correction.

- `JASA_AudibleFieldSweepm_integrated_png.m`  
  Integrated audible-field sweep and compact PNG plotting version.

- `JASA_UltrasonicField.m`  
  Ultrasonic pressure-field calculation.

- `JASA_UltrasonicVelocity.m`  
  Ultrasonic velocity-field calculation.

- `JASA_Plot2DUltral.m`  
  Two-dimensional ultrasonic-field visualization.

## Core Functions

The following local functions are required by the main scripts and are included in this repository:

- `AbsorpAttenCoef.m`  
  Atmospheric absorption coefficient based on ISO 9613-1.

- `solve_kappa0.m`  
  Numerical preparation for the fast Hankel transform.

- `m_FHT.m`  
  m-th order fast Hankel transform implementation.

- `make_source_velocity.m`  
  Source velocity generation for different source profiles and modal orders.

- `calc_ultrasound_field.m`  
  Ultrasonic pressure-field calculation using King/FHT and/or DIM.

- `calc_ultrasound_velocity_field.m`  
  Ultrasonic velocity-field calculation using King/FHT and/or DIM.

- `MyColor.m`  
  Colormap utility used by plotting scripts.

## Requirements

The codes require MATLAB and the local helper functions included in this repository.

The scripts were developed for numerical acoustic-field simulations involving:

- fast Hankel transforms,
- vortex source profiles,
- ultrasonic primary fields,
- nonlinear audible-field generation,
- local-effect correction,
- direct integration comparison,
- convergence and computational-cost analysis.

Parallel computing may be used in some scripts for DIM calculations.

## Typical Workflow

A typical workflow is:

1. Run ultrasonic-field convergence tests:
   ```matlab
   JASA_Calc_convergence_ultra
   ```

2. Run audible-field convergence tests:
   ```matlab
   JASA_Calc_convergence_audio
   ```

3. Run fixed-line King/DIM comparison:
   ```matlab
   JASA_KingDIM_lineCompare
   ```

4. Run computational-cost benchmark:
   ```matlab
   JASA_KingDIM_ComplexityCompare
   ```

5. Run audible-field modal sweep with local-effect correction:
   ```matlab
   JASA_AudibleFieldSweepm
   ```

## Output Files

The scripts may generate result folders containing `.mat`, `.fig`, `.png`, `.pdf`, `.csv`, and `.txt` files. These generated files are not included in this repository by default.

Typical generated folders include:

- `result_*/`
- `data_pu/`
- `AudioTransform_KingLocalEffect_mSweep/`
- `AudioFixed_LocalCorr_*/`
- `JASA_USE/`

These folders are intended for local calculation results and are usually excluded from version control.

## Notes

- The repository is intended mainly for code sharing and reproducible numerical modeling.
- Large cache files and generated result files should not be committed to the repository.
- Some scripts may require substantial memory and computation time, especially DIM-based benchmark calculations.
- The default parameters in each script should be checked before running large-scale simulations.

## License

No license has been specified yet. Add a license file if the repository will be made public or distributed.
