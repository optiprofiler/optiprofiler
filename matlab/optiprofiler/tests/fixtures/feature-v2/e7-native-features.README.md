# e7-native-features.mat

Native `Feature` objects written by the pre-repair candidate e7d5341fb96026ffcd9c782f95081b47e3d4dfab
itself (MATLAB R2026a on syu-ubuntu, 2026-09-14), from values that e7d5341 accepted:

| Variable | Constructed as |
| --- | --- |
| `nan_perturbation` | `Feature('perturbed_x0', 'perturbation_level', NaN)` |
| `negative_perturbation` | `Feature('perturbed_x0', 'perturbation_level', -0.1)` |
| `row_perturbation` | `Feature('perturbed_x0', 'perturbation_level', [0.1 0.2])` |
| `inf_noise` | `Feature('noisy', 'noise_level', Inf)` |
| `nan_mesh` | `Feature('quantized', 'mesh_size', NaN)` |
| `numeric_logical` | `Feature('linearly_transformed', 'rotated', 1)` |
| `int32_digits` | `Feature('truncated', 'significant_digits', int32(3), 'perturbed_trailing_digits', true)` |
| `valid_control` | `Feature('noisy', 'noise_level', 1e-3)` |

Producer: `probe_e7_audit.m` (sha256 4d1b4bb090cd5613df6e46e98dbcd20e1489a5ec5321927df9deedf4c276c32f),
receipt `probe-e7-baseline-attempt2` of the 2026-09-14 debugging round. File sha256 is recorded in that
round's FINAL_REVIEW.md. It is a genuine historical producer fixture; never regenerate it with newer code.
