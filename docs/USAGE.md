# Usage

Open the project folder as MATLAB's Current Folder and run:

```matlab
addpath(genpath(pwd))
```

## Analysis

Run `model/tumor_model_parameters.m` before `analysis/equilibria_analysis.m`. The other scripts in `analysis/` can be run independently.

## Control simulations

For LQR, run `lqr/LQR_design.m` first. Then, for any control method:

1. Run `model/tumor_model_parameters.m`.
2. Set `p.R3 = 4.5` in the workspace.
3. Open the corresponding model and press **Run**.

| Method | Model |
| --- | --- |
| LQR | `lqr/LQR_simulation.slx` |
| Feedback linearization | `feedback-linearization/FBL_simulation.slx` |
| Sliding mode control | `sliding-mode/SMC_simulation.slx` |
| Sliding mode with EKF | `state-estimation/SMC_EKF_simulation.slx` |

## Supporting scripts

- `sliding-mode/lie_derivatives.m`: symbolic control calculations.
- `feedback-linearization/zero_dynamics_phase_portrait.m`: error-dynamics phase portrait.

See the [report](Tumor_growth_model_control.pdf) for the complete study.
