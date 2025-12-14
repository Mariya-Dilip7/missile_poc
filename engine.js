const g0 = 9.80665; // m/s^2

function fmt(x) {
  if (!isFinite(x)) return "–";
  const s = x.toFixed(6);
  return s.replace(/\.?0+$/, "");
}

function readPos(id) {
  const el = document.getElementById(id);
  if (!el) return NaN;
  const v = parseFloat(el.value);
  return v > 0 && isFinite(v) ? v : NaN;
}

function readCommonPos(id) {
  const el = document.getElementById(id);
  if (!el) return NaN;
  const v = parseFloat(el.value);
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
  return mdot * Ve + (pe - pa) * Ae;
}

function computeEngine() {
  const gamma = readCommonPos("commonGamma");
  const R     = readCommonPos("commonR");
  const Tc    = readCommonPos("commonTc");
  const pc    = readCommonPos("commonPc");

  const pe    = readPos("pe");
  const mdot  = readPos("mdot");
  const Ae    = readPos("Ae");
  const paEl  = document.getElementById("pa");
  const paVal = paEl ? parseFloat(paEl.value) : NaN;

  const m0    = readPos("m0");
  const mf    = readPos("mf");

  const errEl = document.getElementById("engineError");
  errEl.textContent = "";

  if ([gamma, R, Tc, pc, pe].some(v => !isFinite(v))) {
    errEl.textContent = "Engine: enter positive common γ, R, T_c, p_c and p_e.";
    return;
  }
  if (pe >= pc) {
    errEl.textContent = "Engine: p_e must be lower than p_c.";
    return;
  }

  const Ve  = calcVe(gamma, R, Tc, pc, pe);
  const CF  = calcCF(gamma, pc, pe);
  const Isp = Ve / g0;

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
    dV = Ve * Math.log(m0 / mf);
  }
  if (isFinite(F) && isFinite(m0)) {
    tw = F / (m0 * g0);
  }
  document.getElementById("dVDisplay").textContent = fmt(dV);
  document.getElementById("twDisplay").textContent = fmt(tw);
}

function resetEngine() {
  const ids = [
    "pe","AeAt",
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
  const inside = (mu * mu) / (sigma * dp);
  if (inside <= 0) return NaN;
  return C * Math.cbrt(inside);
}

function calcAinj(mdot, rho, u) {
  return mdot / (rho * u);
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

  let Ainj  = NaN;
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
  const K = term1 * term2;
  if (K <= 0) return NaN;
  return mdot / (pc * K);
}

function computeNozzle() {
  const gamma = readCommonPos("commonGamma");
  const R     = readCommonPos("commonR");
  const Tc    = readCommonPos("commonTc");
  const pc    = readCommonPos("commonPc");

  const mdot     = readPos("nozMdot");
  const AeAt     = readPos("nozAeAt");
  const thetaDiv = readPos("nozThetaDiv");

  const errEl = document.getElementById("nozError");
  errEl.textContent = "";

  if ([gamma, R, Tc, pc, mdot].some(v => !isFinite(v))) {
    errEl.textContent = "Nozzle: enter common γ, R, T_c, p_c and ṁ.";
    return;
  }

  const At = calcAt(gamma, R, Tc, pc, mdot);
  let Ae   = NaN;
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

// -------- COMMON RESET --------

function resetCommon() {
  ["commonGamma","commonR","commonTc","commonPc"].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.value = "";
  });
}

// -------- RANGE MODE CONFIG (for sliders) --------

const rangeConfig = {
  commonGamma: { min: 1.0, max: 2.0, step: 0.001 },
  commonR:     { min: 100, max: 1000, step: 1 },
  commonTc:    { min: 500, max: 4000, step: 10 },
  commonPc:    { min: 1e5, max: 1e7, step: 1e4 },

  pe:          { min: 1e4, max: 1e6, step: 1000 },
  AeAt:        { min: 1, max: 50, step: 0.1 },
  mdot:        { min: 0.1, max: 50, step: 0.1 },
  Ae:          { min: 0.001, max: 0.5, step: 0.001 },
  pa:          { min: 1e4, max: 1e6, step: 1000 },
  m0:          { min: 1, max: 1000, step: 1 },
  mf:          { min: 0.1, max: 1000, step: 1 },
  tb:          { min: 0.1, max: 200, step: 0.1 },

  injC:        { min: 0.1, max: 2, step: 0.01 },
  injMu:       { min: 1e-5, max: 0.01, step: 1e-5 },
  injSigma:    { min: 0.001, max: 0.2, step: 0.001 },
  injDp:       { min: 1e4, max: 2e6, step: 1e4 },
  injMdot:     { min: 0.1, max: 50, step: 0.1 },
  injRho:      { min: 100, max: 3000, step: 10 },
  injU:        { min: 1, max: 200, step: 1 },
  injHoles:    { min: 1, max: 200, step: 1 },

  nozMdot:     { min: 0.1, max: 50, step: 0.1 },
  nozAeAt:     { min: 1, max: 50, step: 0.1 },
  nozThetaDiv: { min: 1, max: 45, step: 0.5 }
};

// create/hide sliders but keep number inputs

function enableRangeMode(enabled) {
  const configIds = Object.keys(rangeConfig);
  configIds.forEach(id => {
    const numInput = document.getElementById(id);
    if (!numInput) return;

    let slider = document.querySelector(`input[data-slider-for="${id}"]`);

    if (enabled) {
      if (!slider) {
        const cfg = rangeConfig[id];
        slider = document.createElement("input");
        slider.type = "range";
        slider.min = cfg.min;
        slider.max = cfg.max;
        slider.step = cfg.step;
        slider.value = numInput.value !== "" ? numInput.value : cfg.min;
        slider.dataset.sliderFor = id;
        slider.style.width = "100%";
        slider.style.marginTop = "4px";

        numInput.insertAdjacentElement("afterend", slider);

        slider.addEventListener("input", () => {
          numInput.value = slider.value;
          const badge = numInput.nextSibling?.nextSibling;
          if (badge && badge.classList && badge.classList.contains("range-value-badge")) {
            badge.textContent = numInput.value;
          }
        });

        numInput.addEventListener("input", () => {
          if (slider && numInput.value !== "") {
            slider.value = numInput.value;
          }
        });
      }
      slider.style.display = "";
    } else {
      if (slider) slider.style.display = "none";
    }
  });
}

function setInputsEnabled(enabled) {
  const inputs = document.querySelectorAll('input[type="number"], input[data-slider-for]');
  inputs.forEach(input => {
    const id = input.id || input.dataset.sliderFor;
    if (!id || id === "rangeModeToggle" || id === "enableEditToggle") return;
    input.disabled = !enabled;
  });
}

// -------- VALUE BADGES --------

function attachRangeValueBadges() {
  const inputs = document.querySelectorAll('input[type="number"]');
  inputs.forEach(input => {
    const id = input.id;
    if (!id || id === "rangeModeToggle" || id === "enableEditToggle") return;

    if (!input.dataset.hasValueBadge) {
      const badge = document.createElement("span");
      badge.className = "range-value-badge";
      badge.style.marginLeft = "6px";
      badge.style.fontSize = "0.72rem";
      badge.style.color = "#9ca4d1";
      badge.textContent = input.value || "";
      input.insertAdjacentElement("afterend", badge);
      input.dataset.hasValueBadge = "1";

      input.addEventListener("input", () => {
        badge.textContent = input.value;
      });
    }
  });
}

// -------- RESET ALL --------

function resetAll() {
  resetCommon();
  resetEngine();
  resetInjector();
  resetNozzle();

  const rangeToggle = document.getElementById("rangeModeToggle");
  if (rangeToggle) {
    rangeToggle.checked = false;
    enableRangeMode(false);
  }

  const enableToggle = document.getElementById("enableEditToggle");
  if (enableToggle) {
    enableToggle.checked = true;
    setInputsEnabled(true);
  }

  const badges = document.querySelectorAll(".range-value-badge");
  badges.forEach(b => { b.textContent = ""; });
}

// -------- wiring --------

window.addEventListener("DOMContentLoaded", () => {
  document.getElementById("engineComputeBtn").addEventListener("click", computeEngine);
  document.getElementById("injComputeBtn").addEventListener("click", computeInjector);
  document.getElementById("nozComputeBtn").addEventListener("click", computeNozzle);

  document.getElementById("engineResetBtn").addEventListener("click", resetEngine);
  document.getElementById("injResetBtn").addEventListener("click", resetInjector);
  document.getElementById("nozResetBtn").addEventListener("click", resetNozzle);
  document.getElementById("commonResetBtn").addEventListener("click", resetCommon);

  const rangeToggle  = document.getElementById("rangeModeToggle");
  const enableToggle = document.getElementById("enableEditToggle");
  const resetAllBtn  = document.getElementById("resetAllBtn");

  if (rangeToggle) {
    rangeToggle.addEventListener("change", () => {
      enableRangeMode(rangeToggle.checked);
    });
  }

  if (enableToggle) {
    enableToggle.addEventListener("change", () => {
      setInputsEnabled(enableToggle.checked);
    });
  }

  if (resetAllBtn) {
    resetAllBtn.addEventListener("click", resetAll);
  }

  attachRangeValueBadges();
});
