// ---------- ThermodynamicModel: constants ----------
const R = 8.314462618;            // J/mol/K
const MW = 16.043;                // g/mol
const NA = 6.02214076e23;
const TC = 190.564, PC = 45.992, OMEGA = 0.01142;   // methane critical constants (K, bar)
const T_TRIPLE = 90.694;          // K
const LJ_EPS_K = 148.0, LJ_SIG_A = 3.73;            // Lennard-Jones CH4: eps/k (K), sigma (angstrom)
const LJ_EPS_KJ = LJ_EPS_K * R / 1000;              // kJ/mol
const LJ_F_PN = (LJ_EPS_K * 1.380649e-23) / (LJ_SIG_A * 1e-10) * 1e12; // eps/sigma in pN

// ---------- CpModel: RRHO ----------
const C2 = 1.438776877;           // hc/k, cm K
const MODES = [[2917, 1], [1534, 2], [3019, 3], [1306, 3]]; // cm-1, degeneracy (NIST CCCBDB)
function cpRRHO(T) {
  let s = 4 * R;
  for (const [nu, g] of MODES) { const x = C2 * nu / T; if (x < 300) { const e = Math.exp(x); s += R * g * x * x * e / ((e - 1) * (e - 1)); } }
  return s;
}
function hRRHO(T) { // H(T) - H(0 K), J/mol
  let s = 4 * R * T;
  for (const [nu, g] of MODES) { const th = C2 * nu, x = th / T; if (x < 300) s += R * g * th / (Math.exp(x) - 1); }
  return s;
}
// ---------- CpModel: Shomate (NIST WebBook, Chase 1998) ----------
const SH = [
  { lo: 298, hi: 1300, c: [-0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, -74.8731] },
  { lo: 1300, hi: 6000, c: [85.81217, 11.26467, -2.114146, 0.13819, -26.42221, -153.5327, -74.8731] }
];
const SH_LO = 298, SH_HI = 6000;
function shSet(T) { return T <= 1300 ? SH[0].c : SH[1].c; }
function cpShomate(T) { const [A, B, C, D, E] = shSet(T), t = T / 1000; return A + B * t + C * t * t + D * t * t * t + E / (t * t); }
function hShomate(T) { // H(T) - H(298.15), J/mol
  const [A, B, C, D, E, F, H] = shSet(T), t = T / 1000;
  return (A * t + B * t * t / 2 + C * t * t * t / 3 + D * t * t * t * t / 4 - E / t + F - H) * 1000;
}
// ideal-gas enthalpy relative to 298.15 K, valid over the whole slider range
function hIdeal(T) { return T >= 298.15 ? hShomate(T) : hRRHO(T) - hRRHO(298.15); }
function cpIdeal(T) { return T >= 298.15 ? cpShomate(T) : cpRRHO(T); }

// ---------- RealGasModel: Peng-Robinson equation of state ----------
const PR_AC = 0.45724 * R * R * TC * TC / (PC * 1e5), PR_B = 0.07780 * R * TC / (PC * 1e5);
const PR_K = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA * OMEGA;
function cubicRoots(a2, a1, a0) { // z^3 + a2 z^2 + a1 z + a0 = 0, real roots
  const q = (3 * a1 - a2 * a2) / 9, r = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54, d = q * q * q + r * r;
  if (d > 0) { const s = Math.cbrt(r + Math.sqrt(d)), t = Math.cbrt(r - Math.sqrt(d)); return [s + t - a2 / 3]; }
  const th = Math.acos(Math.max(-1, Math.min(1, r / Math.sqrt(-q * q * q)))), m = 2 * Math.sqrt(-q);
  return [m * Math.cos(th / 3) - a2 / 3, m * Math.cos((th + 2 * Math.PI) / 3) - a2 / 3, m * Math.cos((th + 4 * Math.PI) / 3) - a2 / 3];
}
function prState(T, Pbar) {
  const P = Pbar * 1e5, sa = 1 + PR_K * (1 - Math.sqrt(T / TC)), alpha = sa * sa;
  const a = PR_AC * alpha, dadT = -PR_AC * PR_K * sa / Math.sqrt(T * TC);
  const A = a * P / (R * R * T * T), B = PR_B * P / (R * T), s2 = Math.SQRT2;
  const roots = cubicRoots(-(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B * B * B)).filter(z => z > B).sort((x, y) => x - y);
  const lnphi = Z => Z - 1 - Math.log(Z - B) - A / (2 * s2 * B) * Math.log((Z + (1 + s2) * B) / (Z + (1 - s2) * B));
  let Z, phase;
  if (roots.length >= 2) {
    const zl = roots[0], zv = roots[roots.length - 1];
    if (lnphi(zl) < lnphi(zv)) { Z = zl; phase = 'liquid'; } else { Z = zv; phase = 'gas'; }
  } else {
    Z = roots[0];
    const Vc = 0.3074 * R * TC / (PC * 1e5), V = Z * R * T / P;
    phase = T >= TC ? (Pbar >= PC ? 'supercritical' : 'gas') : (V < Vc ? 'liquid' : 'gas');
  }
  if (T < T_TRIPLE) phase = 'solid';
  const L = Math.log((Z + (1 + s2) * B) / (Z + (1 - s2) * B));
  const Hres = R * T * (Z - 1) + (T * dadT - a) / (2 * s2 * PR_B) * L;   // J/mol
  const Ures = Hres - R * T * (Z - 1);
  const Vm = Z * R * T / P;                                             // m3/mol
  return { Z, phase, Hres, Ures, Vm };
}
const isFluidGas = ph => ph === 'gas' || ph === 'supercritical';
function cpResidual(T, P) { // J/mol/K, numerical dHres/dT at constant P
  const h = 0.25, s0 = prState(T, P), a = prState(T - h, P), b = prState(T + h, P);
  if (a.phase === s0.phase && b.phase === s0.phase) return (b.Hres - a.Hres) / (2 * h);
  if (b.phase === s0.phase) return (b.Hres - s0.Hres) / h;
  return (s0.Hres - a.Hres) / h;
}
function cpReal(T, P) { return cpIdeal(T) + cpResidual(T, P); }

// ---------- MolecularModel: kinetic theory + mapping to the simulation ----------
function meanSpeed(T) { return Math.sqrt(8 * R * T / (Math.PI * MW / 1000)); }   // m/s
function meanSpacingA(Vm) { return Math.cbrt(Vm / NA) * 1e10; }                // angstrom
function ljU(r) { const s6 = Math.pow(1 / r, 6); return 4 * (s6 * s6 - s6); } // reduced
function ljF(r) { const s6 = Math.pow(1 / r, 6); return 24 * (2 * s6 * s6 - s6) / r; } // reduced, + = repulsive

// ---------- ThermalLoadModel ----------
function thermalLoad(mdot, T1, T2, P) { // kg/s, K, K, bar -> W
  const n = mdot / (MW / 1000);
  const dShom = (T1 >= SH_LO && T2 >= SH_LO) ? hShomate(T2) - hShomate(T1) : null;
  const dRRHO = hRRHO(T2) - hRRHO(T1);
  const s1 = prState(T1, P), s2 = prState(T2, P);
  const realOK = isFluidGas(s1.phase) && isFluidGas(s2.phase);
  const dReal = realOK ? hIdeal(T2) - hIdeal(T1) + s2.Hres - s1.Hres : null;
  return { n, dShom, dRRHO, dReal, qShom: dShom === null ? null : n * dShom, qRRHO: n * dRRHO, qReal: dReal === null ? null : n * dReal, ph1: s1.phase, ph2: s2.phase };
}
