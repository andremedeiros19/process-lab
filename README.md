# Case 00 · Constant Cp or Cp(T) in a natural-gas preheater

**Question.** In the preliminary duty estimate of a steam methane reformer's natural-gas preheater, up to what outlet temperature is a constant Cp acceptable, and from when must Cp(T) be used?

**Answer.** For a CH₄ preheater fed at 25 °C, a constant Cp taken at 298 K keeps the duty error within 5 % only up to about **110 °C**, whereas Cp evaluated at the **mean temperature** stays within **1.5 %** up to 800 °C. At 30 bar, the ideal-gas assumption itself adds **2.8 %** to the duty (DWSIM, Peng–Robinson) and should be checked in a simulator.

![Error of each constant-Cp shortcut on the heat duty against the outlet temperature](figures/fig1_shortcut_error.png)

## Base case

Pure CH₄, 1 kmol/s, 25 → 370 °C, 30 bar. Reference: Shomate Cp(T) integrated (NIST-JANAF), **15.44 MW**.

| Shortcut | Cp used (kJ/kmol·K) | Duty (MW) | Error |
| --- | ---: | ---: | ---: |
| Cp at the inlet, 298 K | 35.65 | 12.30 | −20.3 % |
| **Cp at the mean temperature, 471 K** | **44.61** | **15.39** | **−0.3 %** |
| Mean of inlet and outlet Cp | 45.16 | 15.58 | +0.9 % |
| Cp at the outlet, 643 K | 54.68 | 18.86 | +22.2 % |

| Real-gas effect (DWSIM) | 1 bar | 10 bar | 20 bar | 30 bar |
| --- | ---: | ---: | ---: | ---: |
| Duty, Peng–Robinson vs ideal gas | +0.09 % | +0.94 % | +1.88 % | **+2.84 %** |

![Real-gas correction against pressure and temperature](figures/fig2_real_gas_dwsim.png)

## Why

Methane's Cp grows from 35.6 to 54.7 J/mol·K between 25 and 370 °C because its bending vibrational modes (θ = 1879 K and 2207 K) switch on inside the exchanger:

$$
C_p(T) = 4R + R\sum_{i=1}^{4} g_i\,\frac{x_i^{2}\,e^{x_i}}{\left(e^{x_i}-1\right)^{2}}, \qquad x_i=\frac{\theta_i}{T}
$$

A Cp at the inlet ignores this. Cp at the mean temperature works because Cp(T) of methane is nearly linear between 300 and 1000 K.

The real-gas excess has a different origin: attraction *between* molecules, which no ideal-gas model includes. It grows with pressure and fades with temperature:

$$
C_p^{\mathrm{real}}(T,P) = C_p^{\mathrm{ig}}(T) + \left(\frac{\partial H^{\mathrm{res}}}{\partial T}\right)_{P}
$$

## Method

| Route | Tool | Agreement |
| --- | --- | --- |
| Statistical mechanics (rigid rotor, harmonic oscillator) | Python, Excel | within 0.1 % of NIST-JANAF |
| Shomate equation (NIST-JANAF) | Excel | reference |
| Real gas, reference equation of state (Setzmann–Wagner) | NIST WebBook, 30 bar | ends at 625 K |
| Ideal gas vs Peng–Robinson | DWSIM | ideal duty within 0.1 % of Shomate; real Cp within 0.4 % of NIST |

Acceptance criterion: |error on the duty| ≤ 5 %, set before the calculations as about one third of a typical 10–15 % exchanger overdesign margin. Tested also at 1 % and 2 %; the decision does not change.

## Contents

| Path | What it is |
| --- | --- |
| [`report.pdf`](report.pdf) | Full report: data, assumptions, equations, results, sensitivity, limitations, references |
| [`figures/`](figures/) | Figures used above |
| [`simulator/`](simulator/) | DWSIM flowsheet: one Heater, ideal gas vs Peng–Robinson |

## Limitations

Pure methane only. Peng–Robinson is less accurate near the critical point and the saturation line. NIST real-gas data stop at 625 K, so the 643 K outlet relies on DWSIM.

## Key references

- Chase, M. W. NIST-JANAF Thermochemical Tables, 4th ed., 1998 (NIST Chemistry WebBook).
- Setzmann, U.; Wagner, W. *J. Phys. Chem. Ref. Data* 20 (1991) 1061–1155. [doi:10.1063/1.555898](https://doi.org/10.1063/1.555898)
- Peng, D.-Y.; Robinson, D. B. *Ind. Eng. Chem. Fundam.* 15 (1976) 59–64.
- Widom, B. *Statistical Mechanics: A Concise Introduction for Chemists*. Cambridge University Press, 2002.

---

André Medeiros · Process Lab · Case 00 · November 2026
