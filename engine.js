const g0 = 9.80665; // m/s^2

function computeVe(gamma, R, Tc, pc, pe) {
  // Ve = sqrt( 2 * gamma / (gamma - 1) * R * Tc * [1 - (pe/pc)^((gamma-1)/gamma)] )
  const pressureRatio = pe / pc;
  if (pressureRatio <= 0 || pressureRatio >= 1) return NaN;

  const exponent = (gamma - 1) / gamma;
  const term = 1 - Math.pow(pressureRatio, exponent);
  if (term <= 0) return NaN;

  const factor = (2 * gamma) / (gamma - 1);
  const inside = factor * R * Tc * term;
  if (inside <= 0) return NaN;

  return Math.sqrt(inside);
}

function computeCF(gamma, pc, pe, AeAt) {
  // Ideal thrust coefficient (no explicit pressure term)
  // CF = sqrt( (2*gamma^2/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1 - (pe/pc)^((gamma-1)/gamma)) )
  const pressureRatio = pe / pc;
  if (pressureRatio <= 0 || pressureRatio >= 1) return NaN;

  const exponent = (gamma - 1) / gamma;
  const bracket = 1 - Math.pow(pressureRatio, exponent);
  if (bracket <= 0) return NaN;

  const base1 = (2 * gamma * gamma) / (gamma - 1);
  const base2 = 2 / (gamma + 1);
  const power2 = (gamma + 1) / (gamma - 1);

  const inside = base1 * Math.pow(base2, power2) * bracket;
  if (inside <= 0) return NaN;

  return Math.sqrt(inside);
}

function computeThrust(mdot, Ve, pe, pa, Ae) {
  // F = mdot * Ve + (pe - pa) * Ae
  return mdot * Ve + (pe - pa) * Ae;
}

// NEW: no rounding/formatting, just raw JS number → string
function formatNumber(x) {
  if (!isFinite(x)) return "–";
  return x.toString(); // preserves full IEEE-754 double precision digits[web:78][web:86]
}

function readPositive(id) {
  const v = parseFloat(document.getElementById(id).value);
  return v > 0 && isFinite(v) ? v : NaN;
}

document.getElementById("computeBtn").addEventListener("click", () => {
  const gamma = readPositive("gamma");
  const R = readPositive("Rexhaust");
  const Tc = readPositive("Tc");
  const pc = readPositive("pc");
  const pe = readPositive("pe");
  const AeAt = readPositive("AeAt");

  const mdotRaw = parseFloat(document.getElementById("mdot").value);
  const AeRaw = parseFloat(document.getElementById("Ae").value);
  const paRaw = parseFloat(document.getElementById("pa").value);

  const errorEl = document.getElementById("error");
  errorEl.textContent = "";

  if ([gamma, R, Tc, pc, pe, AeAt].some((v) => !isFinite(v))) {
    errorEl.textContent =
      "Enter positive values for γ, R, T_c, p_c, p_e and A_e/A_t.";
    return;
  }

  if (pe >= pc) {
    errorEl.textContent = "Exit pressure must be lower than chamber pressure.";
    return;
  }

  const Ve = computeVe(gamma, R, Tc, pc, pe);
  const CF = computeCF(gamma, pc, pe, AeAt);
  const Isp = Ve / g0;

  document.getElementById("VeDisplay").textContent = formatNumber(Ve);
  document.getElementById("IspDisplay").textContent = formatNumber(Isp);
  document.getElementById("CFDisplay").textContent = formatNumber(CF);

  let F = NaN;
  if (isFinite(mdotRaw) && mdotRaw > 0 && isFinite(AeRaw) && AeRaw > 0) {
    const pa = isFinite(paRaw) && paRaw > 0 ? paRaw : 101325;
    F = computeThrust(mdotRaw, Ve, pe, pa, AeRaw);
  }
  document.getElementById("FDisplay").textContent = formatNumber(F);
});

document.getElementById("resetBtn").addEventListener("click", () => {
  document
    .querySelectorAll("input[type='number']")
    .forEach((input) => (input.value = ""));
  document.getElementById("VeDisplay").textContent = "–";
  document.getElementById("IspDisplay").textContent = "–";
  document.getElementById("CFDisplay").textContent = "–";
  document.getElementById("FDisplay").textContent = "–";
  document.getElementById("error").textContent = "";
});
