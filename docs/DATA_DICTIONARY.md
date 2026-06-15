# Data dictionary

## Common time-series columns

| Column | Meaning |
|---|---|
| `t` | nondimensional time |
| `xc` | droplet center position |
| `theta` | angular position, `2*pi*xc/L` |
| `velocity` | droplet velocity |
| `Wbar` | well-mixed vapor concentration |

## Profile columns

| Column | Meaning |
|---|---|
| `x` | periodic interface coordinate |
| `u` | ES concentration |
| `s` | SDS concentration |
| `c` | SDS--ES complex concentration |
| `Wbar` | well-mixed vapor concentration |
| `xc` | droplet center position |

## Perturbation metrics

| Column | Meaning |
|---|---|
| `s0` | SDS reservoir concentration |
| `A_delta_x_late` | late-time displacement activity |
| `V_rms_late` | late-time velocity activity |
