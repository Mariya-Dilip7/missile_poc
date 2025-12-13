const g0 = 9.80665; // m/s^2

function fmt(x) {
  if (!isFinite(x)) return "–";
  const s = x.toFixed(6); // 6 decimals[web:94][web:95]
  return s.replace(/\.?0+$/, "");
}

function readPos(id) {
  const v = parseFloat(document.getElementById(id).value);
  return v > 0 && isFinite(v) ? v : NaN;
}

// -------- ENGINE --------

function calcVe(gamma, R, Tc, pc, pe) {
  const pr = pe / pc;
  if (pr <= 0 || pr >= 1) return NaN;
  const exp = (gamma - 1) / gamma;
  const term = 1 - Math.pow(pr, exp);
  if (term <= 0) return NaN;
  const factor = (2 * gamma) / (gamma - 1);
  const inside = factor * R * Tc * term;
  if (inside <= 0) return NaN;
  return Math.sqrt(inside);
}

function calcCF(gamma, pc, pe) {
  const pr = pe / pc;
  if (pr <= 0 || pr >= 1) return NaN;
  const exp = (gamma - 1) / gamma;
  const bracket = 1 - Math.pow(pr, exp);
  if (bracket <= 0) return NaN;
  const base1 = (2 * gamma * gamma) / (gamma - 1);
  const base2 = 2 / (gamma + 1);
  const power2 = (gamma + 1) / (gamma - 1);
  const inside = base1 * Math.pow(base2, power2) * bracket;
  if (inside <= 0) return NaN;
  return Math.sqrt(inside);
}

function calcThrust(mdot, Ve, pe, pa, Ae) {
  return mdot * Ve + (pe - pa) * Ae; // thrust eq.[web:19][web:21]
}

function computeEngine() {
  const gamma = readPos("gamma");
  const R     = readPos("Rexhaust");
  const Tc    = readPos("Tc");
  const pc    = readPos("pc");
  const pe    = readPos("pe");
  const mdot  = readPos("mdot");
  const Ae    = readPos("Ae");
  const paVal = parseFloat(document.getElementById("pa").value);

  const m0    = readPos("m0");
  const mf    = readPos("mf");

  const errEl = document.getElementById("engineError");
  errEl.textContent = "";

  if ([gamma, R, Tc, pc, pe].some(v => !isFinite(v))) {
    errEl.textContent = "Engine: enter positive γ, R, T_c, p_c, p_e.";
    return;
  }
  if (pe >= pc) {
    errEl.textContent = "Engine: p_e must be lower than p_c.";
    return;
  }

  const Ve = calcVe(gamma, R, Tc, pc, pe);
  const CF = calcCF(gamma, pc, pe);
  const Isp = Ve / g0; // I_sp[web:17][web:35]

  document.getElementById("VeDisplay").textContent  = fmt(Ve);
  document.getElementById("IspDisplay").textContent = fmt(Isp);
  document.getElementById("CFDisplay").textContent  = fmt(CF);

  let F = NaN;
  if (isFinite(mdot) && isFinite(Ae)) {
    const pa = isFinite(paVal) ? paVal : 101325;
    F = calcThrust(mdot, Ve, pe, pa, Ae);
  }
  document.getElementById("FDisplay").textContent = fmt(F);

  let dV = NaN;
  let tw = NaN;
  if (isFinite(m0) && isFinite(mf) && m0 > mf) {
    dV = Ve * Math.log(m0 / mf); // Tsiolkovsky[web:60]
  }
  if (isFinite(F) && isFinite(m0)) {
    tw = F / (m0 * g0);
  }
  document.getElementById("dVDisplay").textContent = fmt(dV);
  document.getElementById("twDisplay").textContent = fmt(tw);
}

function resetEngine() {
  const ids = [
    "gamma","Rexhaust","Tc","pc","pe","AeAt",
    "mdot","Ae","pa","m0","mf","tb"
  ];
  ids.forEach(id => {
    const el = document.getElementById(id);
    if (el) el.value = "";
  });
  document.getElementById("VeDisplay").textContent  = "–";
  document.getElementById("IspDisplay").textContent = "–";
  document.getElementById("CFDisplay").textContent  = "–";
  document.getElementById("FDisplay").textContent   = "–";
  document.getElementById("dVDisplay").textContent  = "–";
  document.getElementById("twDisplay").textContent  = "–";
  document.getElementById("engineError").textContent = "";
}

// -------- INJECTOR --------

function calcD32(C, mu, sigma, dp) {
  const inside = (mu * mu) / (sigma * dp); // SMD correlation base[web:130][web:133]
  if (inside <= 0) return NaN;
  return C * Math.cbrt(inside);
}

function calcAinj(mdot, rho, u) {
  return mdot / (rho * u); // A_inj[web:63]
}

function computeInjector() {
  const C     = readPos("injC");
  const mu    = readPos("injMu");
  const sigma = readPos("injSigma");
  const dp    = readPos("injDp");
  const mdot  = readPos("injMdot");
  const rho   = readPos("injRho");
  const u     = readPos("injU");
  const Nh    = readPos("injHoles");

  const errEl = document.getElementById("injError");
  errEl.textContent = "";

  if ([C, mu, sigma, dp].some(v => !isFinite(v))) {
    errEl.textContent = "Injector: enter C, μ, σ, Δp.";
    return;
  }

  const D32 = calcD32(C, mu, sigma, dp);
  document.getElementById("D32Display").textContent = fmt(D32);

  let Ainj = NaN;
  let Dhole = NaN;
  if ([mdot, rho, u].every(v => isFinite(v))) {
    Ainj = calcAinj(mdot, rho, u);
    if (isFinite(Nh) && Nh > 0) {
      const Ahole = Ainj / Nh;
      Dhole = Math.sqrt((4 * Ahole) / Math.PI);
    }
  }
  document.getElementById("AinjDisplay").textContent  = fmt(Ainj);
  document.getElementById("DholeDisplay").textContent = fmt(Dhole);
}

function resetInjector() {
  const ids = [
    "injC","injMu","injSigma","injDp",
    "injMdot","injRho","injU","injHoles"
  ];
  ids.forEach(id => {
    const el = document.getElementById(id);
    if (el) el.value = "";
  });
  document.getElementById("D32Display").textContent   = "–";
  document.getElementById("AinjDisplay").textContent  = "–";
  document.getElementById("DholeDisplay").textContent = "–";
  document.getElementById("injError").textContent     = "";
}

// -------- NOZZLE --------

function calcAt(gamma, R, Tc, pc, mdot) {
  const term1 = Math.sqrt(gamma / (R * Tc));
  const term2 = Math.pow(2 / (gamma + 1), (gamma + 1) / (2 * (gamma - 1)));
  const K = term1 * term2; // choked-flow factor[web:126][web:132]
  if (K <= 0) return NaN;
  return mdot / (pc * K);
}

function computeNozzle() {
  const gamma = readPos("nozGamma");
  const R     = readPos("nozR");
  const Tc    = readPos("nozTc");
  const pc    = readPos("nozPc");
  const mdot  = readPos("nozMdot");
  const AeAt  = readPos("nozAeAt");
  const thetaDiv = readPos("nozThetaDiv");

  const errEl = document.getElementById("nozError");
  errEl.textContent = "";

  if ([gamma, R, Tc, pc, mdot].some(v => !isFinite(v))) {
    errEl.textContent = "Nozzle: enter γ, R, T_c, p_c, ṁ.";
    return;
  }

  const At = calcAt(gamma, R, Tc, pc, mdot);
  let Ae = NaN;
  if (isFinite(AeAt)) Ae = AeAt * At;

  const Dt = Math.sqrt((4 * At) / Math.PI);
  const De = isFinite(Ae) ? Math.sqrt((4 * Ae) / Math.PI) : NaN;

  let Ldiv = NaN;
  if (isFinite(thetaDiv) && isFinite(De) && isFinite(Dt)) {
    const thetaRad = (thetaDiv * Math.PI) / 180;
    Ldiv = (De - Dt) / (2 * Math.tan(thetaRad));
  }

  document.getElementById("AtDisplay").textContent    = fmt(At);
  document.getElementById("AeNozDisplay").textContent = fmt(Ae);
  document.getElementById("DtDisplay").textContent    = fmt(Dt);
  document.getElementById("DeDisplay").textContent    = fmt(De);
  document.getElementById("LdivDisplay").textContent  = fmt(Ldiv);
}

function resetNozzle() {
  const ids = [
    "nozGamma","nozR","nozTc","nozPc",
    "nozMdot","nozAeAt","nozThetaDiv"
  ];
  ids.forEach(id => {
    const el = document.getElementById(id);
    if (el) el.value = "";
  });
  document.getElementById("AtDisplay").textContent     = "–";
  document.getElementById("AeNozDisplay").textContent  = "–";
  document.getElementById("DtDisplay").textContent     = "–";
  document.getElementById("DeDisplay").textContent     = "–";
  document.getElementById("LdivDisplay").textContent   = "–";
  document.getElementById("nozError").textContent      = "";
}

// -------- wiring --------

window.addEventListener("DOMContentLoaded", () => {
  document.getElementById("engineComputeBtn").addEventListener("click", computeEngine);
  document.getElementById("injComputeBtn").addEventListener("click", computeInjector);
  document.getElementById("nozComputeBtn").addEventListener("click", computeNozzle);

  document.getElementById("engineResetBtn").addEventListener("click", resetEngine);
  document.getElementById("injResetBtn").addEventListener("click", resetInjector);
  document.getElementById("nozResetBtn").addEventListener("click", resetNozzle);
});
