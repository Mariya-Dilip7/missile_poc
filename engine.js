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

function calcCharVel(gamma, R, Tc, pc, mdot, At) {
  // Characteristic velocity: c* = Pc * At / mdot
  if (mdot <= 0 || At <= 0) return NaN;
  return (pc * At) / mdot;
}

function calcAccel(F, m0) {
  // Acceleration: a = F / m0 (or n = a / g0 for g-count)
  if (m0 <= 0 || !isFinite(F)) return NaN;
  return F / m0;
}

function computeEngine() {
  const gamma = readCommonPos("commonGamma");
  const R     = readCommonPos("commonR");
  const Tc    = readCommonPos("commonTc");
  const pc    = readCommonPos("commonPc");

  const pe    = readPos("pe");
  const mdot  = readPos("mdot");
  const Ae    = readPos("Ae");
  const At    = readPos("At");
  const paEl  = document.getElementById("pa");
  const paVal = paEl ? parseFloat(paEl.value) : NaN;

  const m0    = readPos("m0");
  const mf    = readPos("mf");
  const tb    = readPos("tb");
  const m_payload = readPos("m_payload");

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

  // New calculations
  let Cstar = NaN;
  let avgMass = NaN;
  let avgThrust = NaN;
  let accel = NaN;
  let n_accel = NaN;
  let totalImpulse = NaN;
  let mdot_derived = NaN;
  let tb_derived = NaN;
  let net_accel = NaN;
  let apogee = NaN;

  // Characteristic velocity: c* = Pc * At / mdot
  if (isFinite(pc) && isFinite(At) && isFinite(mdot) && mdot > 0) {
    Cstar = (pc * At) / mdot;
  }

  if (isFinite(mdot) && isFinite(tb)) {
    avgMass = m0 - (mdot * tb) / 2;
    avgThrust = isFinite(F) ? F : NaN;
    if (isFinite(avgMass)) {
      accel = avgThrust / avgMass;
      n_accel = accel / g0;
    }
    totalImpulse = F * tb;
  } else if (isFinite(m0) && isFinite(mf) && isFinite(tb)) {
    mdot_derived = (m0 - mf) / tb;
    avgMass = (m0 + mf) / 2;
    avgThrust = isFinite(F) ? F : NaN;
    if (isFinite(avgMass)) {
      accel = avgThrust / avgMass;
      n_accel = accel / g0;
    }
    totalImpulse = F * tb;
  }

  // Net upward acceleration: a_net = a - g0 (during burn)
  if (isFinite(accel)) {
    net_accel = accel - g0;
    
    // Apogee height (simplified): H = (dV^2) / (2 * g0)
    if (isFinite(dV)) {
      apogee = (dV * dV) / (2 * g0);
    }
  }

  document.getElementById("CstarDisplay").textContent = fmt(Cstar);
  document.getElementById("avgMassDisplay").textContent = fmt(avgMass);
  document.getElementById("avgThrustDisplay").textContent = fmt(avgThrust);
  document.getElementById("accelDisplay").textContent = fmt(accel);
  document.getElementById("nAccelDisplay").textContent = fmt(n_accel);
  document.getElementById("totalImpulseDisplay").textContent = fmt(totalImpulse);
  document.getElementById("mdotDisplay").textContent = fmt(mdot_derived);
  document.getElementById("tbDisplay").textContent = fmt(tb_derived);
  document.getElementById("netAccelDisplay").textContent = fmt(net_accel);
  document.getElementById("apogeeDisplay").textContent = fmt(apogee);
}

function resetEngine() {
  const ids = [
    "pe","AeAt",
    "mdot","Ae","At","pa","m0","mf","tb","m_payload"
  ];
  ids.forEach(id => {
    clearInputAndBadge(id);
  });
  document.getElementById("VeDisplay").textContent  = "–";
  document.getElementById("IspDisplay").textContent = "–";
  document.getElementById("CFDisplay").textContent  = "–";
  document.getElementById("FDisplay").textContent   = "–";
  document.getElementById("dVDisplay").textContent  = "–";
  document.getElementById("twDisplay").textContent  = "–";
  document.getElementById("CstarDisplay").textContent = "–";
  document.getElementById("avgMassDisplay").textContent = "–";
  document.getElementById("avgThrustDisplay").textContent = "–";
  document.getElementById("accelDisplay").textContent = "–";
  document.getElementById("nAccelDisplay").textContent = "–";
  document.getElementById("totalImpulseDisplay").textContent = "–";
  document.getElementById("mdotDisplay").textContent = "–";
  document.getElementById("tbDisplay").textContent = "–";
  document.getElementById("netAccelDisplay").textContent = "–";
  document.getElementById("apogeeDisplay").textContent = "–";
  document.getElementById("engineError").textContent = "";
}

// -------- FUEL & OXIDIZER --------

function computeFuelOxidizer() {
  const m0  = readPos("fuelM0");
  const ofRatio = readPos("fuelOFRatio");
  const mdot_fuel = readPos("fuelMdotFuel");
  const mdot_ox = readPos("fuelMdotOx");
  const rho_fuel = readPos("fuelRhoFuel");
  const rho_ox = readPos("fuelRhoOx");
  const A_regression = readPos("fuelARegression");
  const c = readPos("fuelC");
  const n_exp = readPos("fuelN");
  const tb = readPos("fuelTb");
  const pc = readPos("fuelPc");

  const errEl = document.getElementById("fuelError");
  errEl.textContent = "";

  let m_fuel = NaN;
  let m_ox = NaN;
  let V_fuel = NaN;
  let V_ox = NaN;
  let V_tank = NaN;
  let mdot_fuel_derived = NaN;
  let mdot_ox_derived = NaN;
  let regression_rate = NaN;
  let flux = NaN;

  // Calculate fuel and oxidizer masses
  if (isFinite(m0) && isFinite(ofRatio) && ofRatio > 0) {
    m_ox = (m0 * ofRatio) / (1 + ofRatio);
    m_fuel = m0 - m_ox;
  }

  // Calculate volumes
  if (isFinite(m_fuel) && isFinite(rho_fuel) && rho_fuel > 0) {
    V_fuel = m_fuel / rho_fuel;
  }
  if (isFinite(m_ox) && isFinite(rho_ox) && rho_ox > 0) {
    V_ox = m_ox / rho_ox;
    // Tank volume = 1.1 * V_ox (safety margin)
    V_tank = 1.1 * V_ox;
  }

  // Derive mass flow rates from mass and burn time
  if (isFinite(m_fuel) && isFinite(tb) && tb > 0) {
    mdot_fuel_derived = m_fuel / tb;
  }
  if (isFinite(m_ox) && isFinite(tb) && tb > 0) {
    mdot_ox_derived = m_ox / tb;
  }

  // Regression rate: r = c * (Gox)^n where Gox = mdot_ox / A_regression
  if (isFinite(mdot_ox || mdot_ox_derived) && isFinite(A_regression) && 
      isFinite(c) && isFinite(n_exp) && A_regression > 0) {
    const actualMdotOx = isFinite(mdot_ox) ? mdot_ox : mdot_ox_derived;
    flux = actualMdotOx / A_regression;
    regression_rate = c * Math.pow(flux, n_exp);
  }

  document.getElementById("fuelMassDisplay").textContent = fmt(m_fuel);
  document.getElementById("oxMassDisplay").textContent = fmt(m_ox);
  document.getElementById("fuelVolDisplay").textContent = fmt(V_fuel);
  document.getElementById("oxVolDisplay").textContent = fmt(V_ox);
  document.getElementById("tankVolDisplay").textContent = fmt(V_tank);
  document.getElementById("mdotFuelDisplay").textContent = fmt(mdot_fuel_derived);
  document.getElementById("mdotOxDisplay").textContent = fmt(mdot_ox_derived);
  document.getElementById("regressionRateDisplay").textContent = fmt(regression_rate);
  document.getElementById("fluxDisplay").textContent = fmt(flux);
}

function resetFuelOxidizer() {
  const ids = [
    "fuelM0","fuelOFRatio","fuelMdotFuel","fuelMdotOx",
    "fuelRhoFuel","fuelRhoOx","fuelARegression","fuelC","fuelN","fuelTb","fuelPc"
  ];
  ids.forEach(id => {
    clearInputAndBadge(id);
  });
  document.getElementById("fuelMassDisplay").textContent = "–";
  document.getElementById("oxMassDisplay").textContent = "–";
  document.getElementById("fuelVolDisplay").textContent = "–";
  document.getElementById("oxVolDisplay").textContent = "–";
  document.getElementById("tankVolDisplay").textContent = "–";
  document.getElementById("mdotFuelDisplay").textContent = "–";
  document.getElementById("mdotOxDisplay").textContent = "–";
  document.getElementById("regressionRateDisplay").textContent = "–";
  document.getElementById("fluxDisplay").textContent = "–";
  document.getElementById("fuelError").textContent = "";
}

// -------- FUEL GRAIN --------

function computeFuelGrain() {
  const Ro = readPos("grainRo");
  const Ri = readPos("grainRi");
  const H = readPos("grainH");
  const GOX = readPos("grainGOX");
  const c = readPos("grainC");
  const n = readPos("grainN");
  const pc = readPos("grainPc");

  const errEl = document.getElementById("grainError");
  errEl.textContent = "";

  let port_area = NaN;
  let thickness = NaN;
  let regression_rate = NaN;

  // Port area: A = π * Ri^2
  if (isFinite(Ri) && Ri > 0) {
    port_area = Math.PI * Ri * Ri;
  }

  // Grain thickness: t = Ro - Ri
  if (isFinite(Ro) && isFinite(Ri) && Ro > Ri) {
    thickness = Ro - Ri;
  }

  // Regression rate: r = c * (GOX)^n
  if (isFinite(c) && isFinite(GOX) && isFinite(n) && GOX > 0) {
    regression_rate = c * Math.pow(GOX, n);
  }

  document.getElementById("portAreaDisplay").textContent = fmt(port_area);
  document.getElementById("thicknessDisplay").textContent = fmt(thickness);
  document.getElementById("regressionRateGrainDisplay").textContent = fmt(regression_rate);
}

function resetFuelGrain() {
  const ids = ["grainRo","grainRi","grainH","grainGOX","grainC","grainN","grainPc"];
  ids.forEach(id => {
    clearInputAndBadge(id);
  });
  document.getElementById("portAreaDisplay").textContent = "–";
  document.getElementById("thicknessDisplay").textContent = "–";
  document.getElementById("regressionRateGrainDisplay").textContent = "–";
  document.getElementById("grainError").textContent = "";
}

// -------- OXIDIZER TANK --------

function computeOxidizerTank() {
  const r = readPos("tankR");
  const h = readPos("tankH");
  const Cd = readPos("tankCd");
  const A_vent = readPos("tankA");
  const Tr = readPos("tankTr");
  const Po = readPos("tankPo");
  const gamma = readPos("tankGamma");
  const R = readPos("tankR_gas");

  const errEl = document.getElementById("tankError");
  errEl.textContent = "";

  let V_cyl = NaN;
  let V_total = NaN;
  let mdot_relief = NaN;

  // Cylinder volume: V = π * r^2 * h
  if (isFinite(r) && isFinite(h) && r > 0 && h > 0) {
    V_cyl = Math.PI * r * r * h;
  }

  // Total volume (with safety margin): V_total = 1.1 * V_cyl
  if (isFinite(V_cyl)) {
    V_total = 1.1 * V_cyl;
  }

  // Choked relief mass flow rate: mdot = Cd * A * Po * sqrt(γ / (R * Tr)) * ((2/(γ+1))^((γ+1)/(2(γ-1))))
  if (isFinite(Cd) && isFinite(A_vent) && isFinite(Po) && isFinite(gamma) && 
      isFinite(R) && isFinite(Tr) && A_vent > 0 && Tr > 0) {
    const term1 = Math.sqrt(gamma / (R * Tr));
    const exp = (gamma + 1) / (2 * (gamma - 1));
    const term2 = Math.pow(2 / (gamma + 1), exp);
    mdot_relief = Cd * A_vent * Po * term1 * term2;
  }

  document.getElementById("V_cylDisplay").textContent = fmt(V_cyl);
  document.getElementById("V_totalDisplay").textContent = fmt(V_total);
  document.getElementById("mdot_reliefDisplay").textContent = fmt(mdot_relief);
}

function resetOxidizerTank() {
  const ids = ["tankR","tankH","tankCd","tankA","tankTr","tankPo","tankGamma","tankR_gas"];
  ids.forEach(id => {
    clearInputAndBadge(id);
  });
  document.getElementById("V_cylDisplay").textContent = "–";
  document.getElementById("V_totalDisplay").textContent = "–";
  document.getElementById("mdot_reliefDisplay").textContent = "–";
  document.getElementById("tankError").textContent = "";
}

// -------- NOZZLE --------

function calcD32(C, mu, sigma, dp) {
  const inside = (mu * mu) / (sigma * dp);
  if (inside <= 0) return NaN;
  return C * Math.cbrt(inside);
}

function calcAinj(mdot, rho, u) {
  return mdot / (rho * u);
}

function calcSprayAngle(theta_deg, D_spray, L) {
  // Θ = 2 * arctan(D_spray / (2*L))
  if (L <= 0) return NaN;
  const thetaRad = 2 * Math.atan(D_spray / (2 * L));
  return (thetaRad * 180) / Math.PI;
}

function calcIdealVolFlow(dp, rho) {
  // Q_ideal = sqrt(2 * ΔP / ρ)
  if (dp < 0 || rho <= 0) return NaN;
  return Math.sqrt((2 * dp) / rho);
}

function calcDischargeCoeff(mdot_actual, mdot_ideal) {
  // Cd = mdot_actual / mdot_ideal
  if (mdot_ideal <= 0) return NaN;
  return mdot_actual / mdot_ideal;
}

function calcOrificeArea(mdot, Cd, dp, rho) {
  // A = mdot / (Cd * sqrt(2 * ΔP / ρ))
  const denom = Cd * Math.sqrt((2 * dp) / rho);
  if (denom <= 0) return NaN;
  return mdot / denom;
}

function calcOrificeDia(A) {
  // D = sqrt(4*A / π)
  if (A <= 0) return NaN;
  return Math.sqrt((4 * A) / Math.PI);
}

function calcMdotOrifice(mdot, N_holes) {
  // ṁ_orifice = mdot / N_holes
  if (N_holes <= 0) return NaN;
  return mdot / N_holes;
}

function calcInjectionVel(dp, rho) {
  // U_inj = sqrt(2 * ΔP / ρ)
  if (dp < 0 || rho <= 0) return NaN;
  return Math.sqrt((2 * dp) / rho);
}

function calcVolFlowRate(mdot, rho) {
  // Q = mdot / ρ
  if (rho <= 0) return NaN;
  return mdot / rho;
}

function calcSpeedOfSound(gamma, R, T) {
  // a = sqrt(γ * R * T)
  if (T <= 0) return NaN;
  return Math.sqrt(gamma * R * T);
}

function calcMachNumber(velocity, speedOfSound) {
  // M = velocity / a
  if (speedOfSound <= 0) return NaN;
  return velocity / speedOfSound;
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
  
  // New inputs
  const D_spray = readPos("injDSpray");
  const L_spray = readPos("injLSpray");
  const Cd = readPos("injCd");
  const N_holes = readPos("injNHoles");
  const gamma_ox = readPos("injGamma");
  const R_ox = readPos("injROx");
  const T_ox = readPos("injTOx");

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
  let sprayAngle = NaN;
  let Q_ideal = NaN;
  let Cd_calc = NaN;
  let A_orifice = NaN;
  let D_orifice = NaN;
  let mdot_orifice = NaN;
  let U_inj = NaN;
  let Q = NaN;
  let a = NaN;
  let M = NaN;

  if ([mdot, rho, u].every(v => isFinite(v))) {
    Ainj = calcAinj(mdot, rho, u);
    if (isFinite(Nh) && Nh > 0) {
      const Ahole = Ainj / Nh;
      Dhole = Math.sqrt((4 * Ahole) / Math.PI);
    }
    
    // Volumetric flow
    Q = calcVolFlowRate(mdot, rho);
    
    // Ideal volumetric flow from pressure drop
    Q_ideal = calcIdealVolFlow(dp, rho);
    
    // Discharge coefficient if actual mdot known
    if (isFinite(Q_ideal) && isFinite(Q) && Q_ideal > 0) {
      Cd_calc = Q / Q_ideal;
    }
  }

  // Orifice calculations
  if (isFinite(mdot) && isFinite(Cd) && isFinite(dp) && isFinite(rho)) {
    A_orifice = calcOrificeArea(mdot, Cd, dp, rho);
    if (isFinite(A_orifice)) {
      D_orifice = calcOrificeDia(A_orifice);
    }
    
    if (isFinite(N_holes) && N_holes > 0) {
      mdot_orifice = calcMdotOrifice(mdot, N_holes);
    }
  }

  // Injection velocity
  if (isFinite(dp) && isFinite(rho)) {
    U_inj = calcInjectionVel(dp, rho);
  }

  // Spray angle
  if (isFinite(D_spray) && isFinite(L_spray)) {
    sprayAngle = calcSprayAngle(0, D_spray, L_spray);
  }

  // Speed of sound and Mach
  if (isFinite(gamma_ox) && isFinite(R_ox) && isFinite(T_ox)) {
    a = calcSpeedOfSound(gamma_ox, R_ox, T_ox);
    if (isFinite(U_inj) && isFinite(a)) {
      M = calcMachNumber(U_inj, a);
    }
  }

  document.getElementById("AinjDisplay").textContent  = fmt(Ainj);
  document.getElementById("DholeDisplay").textContent = fmt(Dhole);
  document.getElementById("sprayAngleDisplay").textContent = fmt(sprayAngle);
  document.getElementById("Q_idealDisplay").textContent = fmt(Q_ideal);
  document.getElementById("CdDisplay").textContent = fmt(Cd_calc);
  document.getElementById("A_orificeDisplay").textContent = fmt(A_orifice);
  document.getElementById("D_orificeDisplay").textContent = fmt(D_orifice);
  document.getElementById("mdot_orificeDisplay").textContent = fmt(mdot_orifice);
  document.getElementById("U_injDisplay").textContent = fmt(U_inj);
  document.getElementById("QDisplay").textContent = fmt(Q);
  document.getElementById("speedOfSoundDisplay").textContent = fmt(a);
  document.getElementById("MachDisplay").textContent = fmt(M);
}

function resetInjector() {
  const ids = [
    "injC","injMu","injSigma","injDp",
    "injMdot","injRho","injU","injHoles",
    "injDSpray","injLSpray","injCd","injNHoles","injGamma","injROx","injTOx"
  ];
  ids.forEach(id => {
    clearInputAndBadge(id);
  });
  document.getElementById("D32Display").textContent   = "–";
  document.getElementById("AinjDisplay").textContent  = "–";
  document.getElementById("DholeDisplay").textContent = "–";
  document.getElementById("sprayAngleDisplay").textContent = "–";
  document.getElementById("Q_idealDisplay").textContent = "–";
  document.getElementById("CdDisplay").textContent = "–";
  document.getElementById("A_orificeDisplay").textContent = "–";
  document.getElementById("D_orificeDisplay").textContent = "–";
  document.getElementById("mdot_orificeDisplay").textContent = "–";
  document.getElementById("U_injDisplay").textContent = "–";
  document.getElementById("QDisplay").textContent = "–";
  document.getElementById("speedOfSoundDisplay").textContent = "–";
  document.getElementById("MachDisplay").textContent = "–";
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
  const thetaConv = readPos("nozThetaConv");
  const pe = readPos("nozPe");
  const nozzle_type = document.getElementById("nozType")?.value || "conical";

  const errEl = document.getElementById("nozError");
  errEl.textContent = "";

  if ([gamma, R, Tc, pc, mdot].some(v => !isFinite(v))) {
    errEl.textContent = "Nozzle: enter common γ, R, T_c, p_c and ṁ.";
    return;
  }

  const At = calcAt(gamma, R, Tc, pc, mdot);
  let Ae   = NaN;
  let AeAt_calc = NaN;
  
  if (isFinite(AeAt)) {
    Ae = AeAt * At;
    AeAt_calc = AeAt;
  } else if (isFinite(pe)) {
    const pr = pe / pc;
    if (pr > 0 && pr < 1) {
      const exp = (gamma - 1) / gamma;
      const term = Math.pow(2 / (gamma + 1), 1 / exp);
      const bracket = Math.pow(pr, 1 / exp) - Math.pow(pr, (gamma + 1) / (2 * exp));
      const top = Math.sqrt(bracket);
      AeAt_calc = (1 / term) * Math.sqrt(((gamma + 1) / (gamma - 1)) * (2 / (gamma + 1))) * top;
      Ae = AeAt_calc * At;
    }
  }

  const Dt = Math.sqrt((4 * At) / Math.PI);
  const De = isFinite(Ae) ? Math.sqrt((4 * Ae) / Math.PI) : NaN;

  let Ldiv = NaN;
  let Ldiv_bell = NaN;
  
  if (isFinite(thetaDiv) && isFinite(De) && isFinite(Dt) && isFinite(AeAt_calc)) {
    const thetaRad = (thetaDiv * Math.PI) / 180;
    
    // Conical nozzle (empirical rule): Ldiv = (Dt / 2) * sqrt(AeAt - 1) / tan(θ)
    if (nozzle_type === "conical") {
      Ldiv = (Dt / 2) * Math.sqrt(AeAt_calc - 1) / Math.tan(thetaRad);
    }
    
    // Bell nozzle (Roa approximation): more complex, using geometric method
    if (nozzle_type === "bell") {
      // Simplified: Ldiv = 1.5 * (Dt / 2) * sqrt(AeAt - 1) / tan(θ)
      Ldiv_bell = 1.5 * (Dt / 2) * Math.sqrt(AeAt_calc - 1) / Math.tan(thetaRad);
    }
  }

  let Lconv = NaN;
  if (isFinite(thetaConv) && isFinite(Dt)) {
    const thetaRad = (thetaConv * Math.PI) / 180;
    // Convergence length (empirical): Lconv = (Dt / 2) / tan(θ)
    Lconv = (Dt / 2) / Math.tan(thetaRad);
  }

  let theta_conv_half = NaN;
  if (isFinite(Lconv) && isFinite(Dt)) {
    // θ_conv = arctan(Dt / (2 * Lconv))
    theta_conv_half = Math.atan((Dt / 2) / Lconv) * (180 / Math.PI);
  }

  document.getElementById("AtDisplay").textContent    = fmt(At);
  document.getElementById("AeNozDisplay").textContent = fmt(Ae);
  document.getElementById("AeAtDisplay").textContent  = fmt(AeAt_calc);
  document.getElementById("DtDisplay").textContent    = fmt(Dt);
  document.getElementById("DeDisplay").textContent    = fmt(De);
  document.getElementById("LdivDisplay").textContent  = fmt(nozzle_type === "conical" ? Ldiv : Ldiv_bell);
  document.getElementById("LconvDisplay").textContent = fmt(Lconv);
  document.getElementById("theta_convDisplay").textContent = fmt(theta_conv_half);
}

function resetNozzle() {
  const ids = [
    "nozMdot","nozAeAt","nozThetaDiv","nozThetaConv","nozPe","nozType"
  ];
  ids.forEach(id => {
    clearInputAndBadge(id);
  });
  document.getElementById("AtDisplay").textContent     = "–";
  document.getElementById("AeNozDisplay").textContent  = "–";
  document.getElementById("AeAtDisplay").textContent   = "–";
  document.getElementById("DtDisplay").textContent     = "–";
  document.getElementById("DeDisplay").textContent     = "–";
  document.getElementById("LdivDisplay").textContent   = "–";
  document.getElementById("LconvDisplay").textContent  = "–";
  document.getElementById("theta_convDisplay").textContent  = "–";
  document.getElementById("nozError").textContent      = "";
}

// -------- COMMON RESET --------

// -------- RESET FUNCTIONS --------

function clearInputAndBadge(id) {
  const el = document.getElementById(id);
  if (el) {
    el.value = "";
    
    // Clear the range input if it exists
    const rangeSlider = document.querySelector(`input[data-slider-for="${id}"]`);
    if (rangeSlider) {
      rangeSlider.value = "";
    }
    
    // Clear the badge - it could be after the number input or after the range slider
    let badge = el.nextElementSibling;
    while (badge) {
      if (badge.classList && badge.classList.contains("range-value-badge")) {
        badge.textContent = "";
        break;
      }
      badge = badge.nextElementSibling;
    }
  }
}

function resetCommon() {
  ["commonGamma","commonR","commonTc","commonPc"].forEach(id => {
    clearInputAndBadge(id);
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
  At:          { min: 0.0001, max: 0.1, step: 0.0001 },
  pa:          { min: 1e4, max: 1e6, step: 1000 },
  m0:          { min: 1, max: 1000, step: 1 },
  mf:          { min: 0.1, max: 1000, step: 1 },
  tb:          { min: 0.1, max: 200, step: 0.1 },
  m_payload:   { min: 1, max: 500, step: 1 },

  injC:        { min: 0.1, max: 2, step: 0.01 },
  injMu:       { min: 1e-5, max: 0.01, step: 1e-5 },
  injSigma:    { min: 0.001, max: 0.2, step: 0.001 },
  injDp:       { min: 1e4, max: 2e6, step: 1e4 },
  injMdot:     { min: 0.1, max: 50, step: 0.1 },
  injRho:      { min: 100, max: 3000, step: 10 },
  injU:        { min: 1, max: 200, step: 1 },
  injHoles:    { min: 1, max: 200, step: 1 },
  injDSpray:   { min: 0.001, max: 1, step: 0.001 },
  injLSpray:   { min: 0.01, max: 1, step: 0.01 },
  injCd:       { min: 0.5, max: 1, step: 0.01 },
  injNHoles:   { min: 1, max: 500, step: 1 },
  injGamma:    { min: 1.0, max: 1.5, step: 0.01 },
  injROx:      { min: 100, max: 500, step: 10 },
  injTOx:      { min: 200, max: 400, step: 5 },

  nozMdot:     { min: 0.1, max: 50, step: 0.1 },
  nozAeAt:     { min: 1, max: 50, step: 0.1 },
  nozThetaDiv: { min: 1, max: 45, step: 0.5 },
  nozThetaConv: { min: 1, max: 90, step: 0.5 },
  nozPe:       { min: 1e4, max: 1e6, step: 1000 },

  fuelM0:      { min: 1, max: 1000, step: 1 },
  fuelOFRatio: { min: 0.1, max: 10, step: 0.1 },
  fuelMdotFuel: { min: 0.1, max: 50, step: 0.1 },
  fuelMdotOx:   { min: 0.1, max: 50, step: 0.1 },
  fuelRhoFuel:  { min: 100, max: 2000, step: 10 },
  fuelRhoOx:    { min: 100, max: 2000, step: 10 },
  fuelARegression: { min: 0.001, max: 1, step: 0.001 },
  fuelC:        { min: 0.0001, max: 0.01, step: 0.0001 },
  fuelN:        { min: 0.5, max: 1.5, step: 0.05 },
  fuelTb:       { min: 0.1, max: 200, step: 0.1 },
  fuelPc:       { min: 1e5, max: 1e7, step: 1e4 },

  grainRo:      { min: 0.01, max: 0.5, step: 0.001 },
  grainRi:      { min: 0.001, max: 0.4, step: 0.001 },
  grainH:       { min: 0.01, max: 2, step: 0.01 },
  grainGOX:     { min: 1, max: 1000, step: 10 },
  grainC:       { min: 0.0001, max: 0.01, step: 0.0001 },
  grainN:       { min: 0.5, max: 1.5, step: 0.05 },
  grainPc:      { min: 1e5, max: 1e7, step: 1e4 },

  tankR:        { min: 0.01, max: 1, step: 0.01 },
  tankH:        { min: 0.1, max: 2, step: 0.1 },
  tankCd:       { min: 0.5, max: 1, step: 0.01 },
  tankA:        { min: 0.0001, max: 0.1, step: 0.0001 },
  tankTr:       { min: 200, max: 400, step: 5 },
  tankPo:       { min: 1e5, max: 1e7, step: 1e4 },
  tankGamma:    { min: 1.0, max: 1.5, step: 0.01 },
  tankR_gas:    { min: 100, max: 500, step: 10 }
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

// -------- LOAD EXAMPLE VALUES --------

function loadExampleValues() {
  const examples = {
    // Common parameters
    commonGamma: 1.22,
    commonR: 285.6,
    commonTc: 3000,
    commonPc: 4000000,

    // Engine
    pe: 101325,
    AeAt: 10,
    mdot: 5,
    Ae: 0.0403,
    At: 0.004,
    pa: 101325,
    m0: 200,
    mf: 100,
    tb: 20,
    m_payload: 50,

    // Injector
    injC: 0.5,
    injMu: 0.0015,
    injSigma: 0.02,
    injDp: 500000,
    injMdot: 4.3,
    injRho: 1140,
    injU: 20,
    injHoles: 16,
    injDSpray: 0.1,
    injLSpray: 0.5,
    injCd: 0.8,
    injNHoles: 32,
    injGamma: 1.4,
    injROx: 300,
    injTOx: 300,

    // Nozzle
    nozMdot: 5,
    nozAeAt: 10,
    nozPe: 101325,
    nozThetaDiv: 15,
    nozThetaConv: 45,

    // Fuel & Oxidizer
    fuelM0: 200,
    fuelOFRatio: 2.0,
    fuelMdotFuel: 1.5,
    fuelMdotOx: 3.0,
    fuelRhoFuel: 880,
    fuelRhoOx: 1400,
    fuelARegression: 0.05,
    fuelC: 0.00334,
    fuelN: 0.65,
    fuelTb: 20,
    fuelPc: 4400000,

    // Fuel Grain
    grainRo: 0.15,
    grainRi: 0.05,
    grainH: 1.0,
    grainGOX: 100,
    grainC: 0.0015,
    grainN: 0.65,
    grainPc: 4400000,

    // Oxidizer Tank
    tankR: 0.2,
    tankH: 1.0,
    tankCd: 0.8,
    tankA: 0.001,
    tankTr: 300,
    tankPo: 4400000,
    tankGamma: 1.4,
    tankR_gas: 300
  };

  // Set all example values
  for (const [id, value] of Object.entries(examples)) {
    const el = document.getElementById(id);
    if (el) {
      el.value = value;
      el.dispatchEvent(new Event("input", { bubbles: true }));
    }
  }
}

// -------- CALCULATE ALL --------

function calculateAll() {
  // Calculate all sections in sequence
  computeEngine();
  setTimeout(() => {
    computeInjector();
  }, 100);
  
  setTimeout(() => {
    computeNozzle();
  }, 200);
  
  setTimeout(() => {
    computeFuelOxidizer();
  }, 300);
  
  setTimeout(() => {
    computeFuelGrain();
  }, 400);
  
  setTimeout(() => {
    computeOxidizerTank();
  }, 500);
}

// -------- RESET ALL --------

function resetAll() {
  resetCommon();
  resetEngine();
  resetInjector();
  resetNozzle();
  resetFuelOxidizer();
  resetFuelGrain();
  resetOxidizerTank();

  // Reset all range sliders to minimum position
  const rangeSliders = document.querySelectorAll('input[data-slider-for]');
  rangeSliders.forEach(slider => {
    const id = slider.dataset.sliderFor;
    if (id && rangeConfig[id]) {
      slider.value = rangeConfig[id].min;
    }
  });

  // Reset all result displays
  const resultValues = document.querySelectorAll(".result-value");
  resultValues.forEach(rv => { rv.textContent = "–"; });

  // Reset all error messages
  const errors = document.querySelectorAll(".error");
  errors.forEach(e => { e.textContent = ""; });

  // Reset checkboxes to default state - this will hide the range sliders
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
}

// -------- wiring --------

window.addEventListener("DOMContentLoaded", () => {
  // Add collapsible section functionality
  const sectionTitles = document.querySelectorAll(".section-title");
  sectionTitles.forEach(title => {
    title.addEventListener("click", () => {
      const content = title.nextElementSibling;
      const toggle = title.querySelector(".section-toggle");
      if (content) {
        content.classList.toggle("collapsed");
        toggle.textContent = content.classList.contains("collapsed") ? "+" : "−";
      }
    });
  });

  document.getElementById("engineComputeBtn").addEventListener("click", computeEngine);
  document.getElementById("injComputeBtn").addEventListener("click", computeInjector);
  document.getElementById("nozComputeBtn").addEventListener("click", computeNozzle);
  document.getElementById("fuelComputeBtn").addEventListener("click", computeFuelOxidizer);
  document.getElementById("grainComputeBtn").addEventListener("click", computeFuelGrain);
  document.getElementById("tankComputeBtn").addEventListener("click", computeOxidizerTank);

  document.getElementById("engineResetBtn").addEventListener("click", resetEngine);
  document.getElementById("injResetBtn").addEventListener("click", resetInjector);
  document.getElementById("nozResetBtn").addEventListener("click", resetNozzle);
  document.getElementById("fuelResetBtn").addEventListener("click", resetFuelOxidizer);
  document.getElementById("grainResetBtn").addEventListener("click", resetFuelGrain);
  document.getElementById("tankResetBtn").addEventListener("click", resetOxidizerTank);
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

  const loadExampleBtn = document.getElementById("loadExampleBtn");
  if (loadExampleBtn) {
    loadExampleBtn.addEventListener("click", loadExampleValues);
  }

  const calculateAllBtn = document.getElementById("calculateAllBtn");
  if (calculateAllBtn) {
    calculateAllBtn.addEventListener("click", () => {
      calculateAll();
      setTimeout(() => {
        updateVisualizations();
      }, 600);
    });
  }

  // -------- VISUALIZATIONS --------

  let charts = {};

  function updateVisualizations() {
    updateEnginePerformanceChart();
    updateThrustAccelChart();
    updateMassVelChart();
    updatePropellantChart();
    updateNozzleDiagram();
    updateTankCapacityChart();
  }

  function updateEnginePerformanceChart() {
    const Ve = parseFloat(document.getElementById("VeDisplay")?.textContent || 0);
    const Isp = parseFloat(document.getElementById("IspDisplay")?.textContent || 0);
    const CF = parseFloat(document.getElementById("CFDisplay")?.textContent || 0);
    const F = parseFloat(document.getElementById("FDisplay")?.textContent || 0);

    const ctx = document.getElementById("enginePerformanceChart");
    if (!ctx) return;

    if (charts.enginePerf) charts.enginePerf.destroy();

    charts.enginePerf = new Chart(ctx, {
      type: "bar",
      data: {
        labels: ["Ve (m/s)", "Isp (s)", "CF", "F (kN)"],
        datasets: [{
          label: "Engine Metrics",
          data: [Ve, Isp, CF, F / 1000],
          backgroundColor: [
            "rgba(100, 217, 255, 0.8)",
            "rgba(91, 136, 255, 0.8)",
            "rgba(100, 217, 255, 0.6)",
            "rgba(91, 136, 255, 0.6)"
          ],
          borderColor: "#64d9ff",
          borderWidth: 2,
          borderRadius: 8
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: true,
        indexAxis: "y",
        plugins: {
          legend: { display: false },
          tooltip: {
            backgroundColor: "rgba(15, 20, 40, 0.9)",
            titleColor: "#64d9ff",
            bodyColor: "#e8eef7",
            borderColor: "#64d9ff",
            borderWidth: 1
          }
        },
        scales: {
          x: {
            ticks: { color: "#b0b8d8" },
            grid: { color: "rgba(100, 217, 255, 0.1)" }
          },
          y: {
            ticks: { color: "#b0b8d8" }
          }
        }
      }
    });
  }

  function updateThrustAccelChart() {
    const F = parseFloat(document.getElementById("FDisplay")?.textContent || 0);
    const accel = parseFloat(document.getElementById("accelDisplay")?.textContent || 0);
    const nAccel = parseFloat(document.getElementById("nAccelDisplay")?.textContent || 0);

    const ctx = document.getElementById("thrustAccelChart");
    if (!ctx) return;

    if (charts.thrustAccel) charts.thrustAccel.destroy();

    // Simulate thrust curve over burn time
    const burnTime = parseFloat(document.getElementById("tb")?.value || 20);
    const timePoints = [];
    const thrustPoints = [];
    const accelPoints = [];

    for (let t = 0; t <= burnTime; t += burnTime / 20) {
      timePoints.push(t.toFixed(1));
      thrustPoints.push((F * (1 - 0.3 * Math.sin(Math.PI * t / burnTime))).toFixed(0));
      accelPoints.push((accel * (1 - 0.2 * Math.sin(Math.PI * t / burnTime))).toFixed(2));
    }

    charts.thrustAccel = new Chart(ctx, {
      type: "line",
      data: {
        labels: timePoints,
        datasets: [
          {
            label: "Thrust (N)",
            data: thrustPoints,
            borderColor: "#64d9ff",
            backgroundColor: "rgba(100, 217, 255, 0.1)",
            tension: 0.3,
            yAxisID: "y",
            borderWidth: 2
          },
          {
            label: "Acceleration (m/s²)",
            data: accelPoints,
            borderColor: "#ff9db8",
            backgroundColor: "rgba(255, 157, 184, 0.1)",
            tension: 0.3,
            yAxisID: "y1",
            borderWidth: 2
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: true,
        interaction: { mode: "index", intersect: false },
        plugins: {
          legend: {
            labels: { color: "#b0b8d8", boxWidth: 12 },
            position: "top"
          },
          tooltip: {
            backgroundColor: "rgba(15, 20, 40, 0.9)",
            titleColor: "#64d9ff",
            bodyColor: "#e8eef7",
            borderColor: "#64d9ff",
            borderWidth: 1
          }
        },
        scales: {
          x: {
            ticks: { color: "#b0b8d8" },
            grid: { color: "rgba(100, 217, 255, 0.1)" }
          },
          y: {
            ticks: { color: "#b0b8d8" },
            grid: { color: "rgba(100, 217, 255, 0.1)" }
          },
          y1: {
            type: "linear",
            position: "right",
            ticks: { color: "#ff9db8" },
            grid: { drawOnChartArea: false }
          }
        }
      }
    });
  }

  function updateMassVelChart() {
    const m0 = parseFloat(document.getElementById("m0")?.value || 0);
    const mf = parseFloat(document.getElementById("mf")?.value || 0);
    const dv = parseFloat(document.getElementById("dVDisplay")?.textContent || 0);

    const ctx = document.getElementById("massVelChart");
    if (!ctx) return;

    if (charts.massVel) charts.massVel.destroy();

    const burnTime = parseFloat(document.getElementById("tb")?.value || 20);
    const timePoints = [];
    const massPoints = [];
    const velPoints = [];

    for (let t = 0; t <= burnTime; t += burnTime / 20) {
      const ratio = t / burnTime;
      timePoints.push(t.toFixed(1));
      massPoints.push((m0 - (m0 - mf) * ratio).toFixed(1));
      velPoints.push((dv * ratio).toFixed(0));
    }

    charts.massVel = new Chart(ctx, {
      type: "line",
      data: {
        labels: timePoints,
        datasets: [
          {
            label: "Mass (kg)",
            data: massPoints,
            borderColor: "#a8d8ff",
            backgroundColor: "rgba(168, 216, 255, 0.1)",
            tension: 0.3,
            yAxisID: "y",
            borderWidth: 2
          },
          {
            label: "Velocity (m/s)",
            data: velPoints,
            borderColor: "#ffc857",
            backgroundColor: "rgba(255, 200, 87, 0.1)",
            tension: 0.3,
            yAxisID: "y1",
            borderWidth: 2
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: true,
        interaction: { mode: "index", intersect: false },
        plugins: {
          legend: {
            labels: { color: "#b0b8d8", boxWidth: 12 },
            position: "top"
          },
          tooltip: {
            backgroundColor: "rgba(15, 20, 40, 0.9)",
            titleColor: "#64d9ff",
            bodyColor: "#e8eef7",
            borderColor: "#64d9ff",
            borderWidth: 1
          }
        },
        scales: {
          x: {
            ticks: { color: "#b0b8d8" },
            grid: { color: "rgba(100, 217, 255, 0.1)" }
          },
          y: {
            ticks: { color: "#b0b8d8" },
            grid: { color: "rgba(100, 217, 255, 0.1)" }
          },
          y1: {
            type: "linear",
            position: "right",
            ticks: { color: "#ffc857" },
            grid: { drawOnChartArea: false }
          }
        }
      }
    });
  }

  function updatePropellantChart() {
    const fuelMass = parseFloat(document.getElementById("fuelMassDisplay")?.textContent || 0);
    const oxMass = parseFloat(document.getElementById("oxMassDisplay")?.textContent || 0);
    const payloadMass = parseFloat(document.getElementById("m_payload")?.value || 0);
    const dryMass = parseFloat(document.getElementById("mf")?.value || 0) - fuelMass - oxMass;

    const ctx = document.getElementById("propellantChart");
    if (!ctx) return;

    if (charts.propellant) charts.propellant.destroy();

    charts.propellant = new Chart(ctx, {
      type: "doughnut",
      data: {
        labels: ["Fuel", "Oxidizer", "Payload", "Dry Mass"],
        datasets: [{
          data: [fuelMass, oxMass, payloadMass, dryMass],
          backgroundColor: [
            "rgba(255, 200, 87, 0.8)",
            "rgba(100, 217, 255, 0.8)",
            "rgba(255, 157, 184, 0.8)",
            "rgba(150, 170, 200, 0.8)"
          ],
          borderColor: "#1a1f2e",
          borderWidth: 2
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: true,
        plugins: {
          legend: {
            labels: { color: "#b0b8d8", boxWidth: 12, padding: 12 },
            position: "bottom"
          },
          tooltip: {
            backgroundColor: "rgba(15, 20, 40, 0.9)",
            titleColor: "#64d9ff",
            bodyColor: "#e8eef7",
            borderColor: "#64d9ff",
            borderWidth: 1,
            callbacks: {
              label: function(context) {
                return context.label + ": " + context.parsed + " kg";
              }
            }
          }
        }
      }
    });
  }

  function updateNozzleDiagram() {
    const svg = document.getElementById("nozzleDiagram");
    if (!svg) return;

    const Dt = parseFloat(document.getElementById("DtDisplay")?.textContent || 0.05);
    const De = parseFloat(document.getElementById("DeDisplay")?.textContent || 0.15);
    const Ldiv = parseFloat(document.getElementById("LdivDisplay")?.textContent || 0.3);
    const Lconv = parseFloat(document.getElementById("LconvDisplay")?.textContent || 0.2);

    // Clear previous content
    svg.innerHTML = "";

    const width = 600, height = 300;
    const scale = 300;
    const centerY = height / 2;

    // Throat position with convergent section
    const chamberX = 40;
    const throatX = chamberX + Lconv * scale;
    const exitX = throatX + Ldiv * scale;

    // Draw chamber (larger)
    const chamberRect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    chamberRect.setAttribute("x", chamberX.toString());
    chamberRect.setAttribute("y", (centerY - Dt * scale / 2).toString());
    chamberRect.setAttribute("width", (Lconv * scale).toString());
    chamberRect.setAttribute("height", (Dt * scale).toString());
    chamberRect.setAttribute("fill", "rgba(100, 217, 255, 0.3)");
    chamberRect.setAttribute("stroke", "#64d9ff");
    chamberRect.setAttribute("stroke-width", "2.5");
    svg.appendChild(chamberRect);

    // Draw throat
    const throatRect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    throatRect.setAttribute("x", throatX.toString());
    throatRect.setAttribute("y", (centerY - Dt * scale / 2).toString());
    throatRect.setAttribute("width", "8");
    throatRect.setAttribute("height", (Dt * scale).toString());
    throatRect.setAttribute("fill", "rgba(255, 200, 87, 0.8)");
    throatRect.setAttribute("stroke", "#ffc857");
    throatRect.setAttribute("stroke-width", "2.5");
    svg.appendChild(throatRect);

    // Draw divergent nozzle (expansion section)
    const divergentPath = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
    const topY = centerY - De * scale / 2;
    const botY = centerY + De * scale / 2;
    const throatTopY = centerY - Dt * scale / 2;
    const throatBotY = centerY + Dt * scale / 2;
    
    divergentPath.setAttribute("points", `${throatX + 8},${throatTopY} ${exitX},${topY} ${exitX},${botY} ${throatX + 8},${throatBotY}`);
    divergentPath.setAttribute("fill", "rgba(100, 217, 255, 0.2)");
    divergentPath.setAttribute("stroke", "#64d9ff");
    divergentPath.setAttribute("stroke-width", "2.5");
    svg.appendChild(divergentPath);

    // Draw centerline
    const centerLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
    centerLine.setAttribute("x1", "20");
    centerLine.setAttribute("y1", centerY.toString());
    centerLine.setAttribute("x2", exitX.toString());
    centerLine.setAttribute("y2", centerY.toString());
    centerLine.setAttribute("stroke", "rgba(100, 217, 255, 0.3)");
    centerLine.setAttribute("stroke-width", "1");
    centerLine.setAttribute("stroke-dasharray", "4,4");
    svg.appendChild(centerLine);

    // Add dimension lines and labels
    // Throat diameter label
    const throatDimLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
    throatDimLine.setAttribute("x1", (throatX - 15).toString());
    throatDimLine.setAttribute("y1", throatTopY.toString());
    throatDimLine.setAttribute("x2", (throatX - 15).toString());
    throatDimLine.setAttribute("y2", throatBotY.toString());
    throatDimLine.setAttribute("stroke", "#ffc857");
    throatDimLine.setAttribute("stroke-width", "1.5");
    svg.appendChild(throatDimLine);

    // Exit diameter label
    const exitDimLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
    exitDimLine.setAttribute("x1", (exitX + 15).toString());
    exitDimLine.setAttribute("y1", topY.toString());
    exitDimLine.setAttribute("x2", (exitX + 15).toString());
    exitDimLine.setAttribute("y2", botY.toString());
    exitDimLine.setAttribute("stroke", "#64d9ff");
    exitDimLine.setAttribute("stroke-width", "1.5");
    svg.appendChild(exitDimLine);

    // Section labels with background
    function createLabel(x, y, lines, color, bgColor) {
      // Handle both single string and array of strings
      const textLines = Array.isArray(lines) ? lines : [lines];
      const bgHeight = textLines.length > 1 ? 28 : 22;
      
      // Background rectangle
      const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      bg.setAttribute("x", (x - 45).toString());
      bg.setAttribute("y", (y - bgHeight/2).toString());
      bg.setAttribute("width", "90");
      bg.setAttribute("height", bgHeight.toString());
      bg.setAttribute("fill", bgColor);
      bg.setAttribute("stroke", color);
      bg.setAttribute("stroke-width", "1.5");
      bg.setAttribute("rx", "4");
      svg.appendChild(bg);

      // Create text group for multi-line support
      const textGroup = document.createElementNS("http://www.w3.org/2000/svg", "text");
      textGroup.setAttribute("x", x.toString());
      textGroup.setAttribute("y", y.toString());
      textGroup.setAttribute("fill", color);
      textGroup.setAttribute("font-size", "12");
      textGroup.setAttribute("font-weight", "700");
      textGroup.setAttribute("text-anchor", "middle");
      textGroup.setAttribute("dominant-baseline", "middle");
      
      textLines.forEach((line, index) => {
        const tspan = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
        tspan.setAttribute("x", x.toString());
        tspan.setAttribute("dy", index === 0 ? "0" : "14");
        tspan.textContent = line;
        textGroup.appendChild(tspan);
      });
      
      svg.appendChild(textGroup);
    }

    // Chamber label
    createLabel(chamberX + Lconv * scale / 2, 30, "CHAMBER", "#64d9ff", "rgba(100, 217, 255, 0.1)");

    // Throat label
    createLabel(throatX + 4, centerY + 50, [`THROAT`, `Dt=${Dt.toFixed(3)}m`], "#ffc857", "rgba(255, 200, 87, 0.1)");

    // Nozzle label
    createLabel(throatX + Ldiv * scale / 2 + 4, centerY - 60, [`DIVERGENT`, `NOZZLE`, `De=${De.toFixed(3)}m`], "#64d9ff", "rgba(100, 217, 255, 0.1)");

    // Expansion ratio label
    const expansionRatio = (De / Dt).toFixed(2);
    const ratioLabel = document.createElementNS("http://www.w3.org/2000/svg", "text");
    ratioLabel.setAttribute("x", "300");
    ratioLabel.setAttribute("y", "280");
    ratioLabel.setAttribute("fill", "#a8d8ff");
    ratioLabel.setAttribute("font-size", "12");
    ratioLabel.setAttribute("text-anchor", "middle");
    ratioLabel.setAttribute("font-weight", "600");
    ratioLabel.textContent = `Expansion Ratio = ${expansionRatio}`;
    svg.appendChild(ratioLabel);

    // Add flow direction arrow
    const arrowGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
    const arrowPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
    arrowPath.setAttribute("d", `M ${exitX - 30} 20 L ${exitX} 20 L ${exitX - 8} 14 M ${exitX} 20 L ${exitX - 8} 26`);
    arrowPath.setAttribute("stroke", "#7af0c2");
    arrowPath.setAttribute("stroke-width", "2");
    arrowPath.setAttribute("fill", "none");
    arrowPath.setAttribute("stroke-linecap", "round");
    arrowPath.setAttribute("stroke-linejoin", "round");
    arrowGroup.appendChild(arrowPath);

    const flowLabel = document.createElementNS("http://www.w3.org/2000/svg", "text");
    flowLabel.setAttribute("x", (exitX - 20).toString());
    flowLabel.setAttribute("y", "12");
    flowLabel.setAttribute("fill", "#7af0c2");
    flowLabel.setAttribute("font-size", "11");
    flowLabel.setAttribute("font-weight", "600");
    flowLabel.textContent = "FLOW";
    arrowGroup.appendChild(flowLabel);

    svg.appendChild(arrowGroup);
  }

  function updateTankCapacityChart() {
    const V_cyl = parseFloat(document.getElementById("V_cylDisplay")?.textContent || 0);
    const V_total = parseFloat(document.getElementById("V_totalDisplay")?.textContent || 0);
    const oxVol = parseFloat(document.getElementById("oxVolDisplay")?.textContent || 0);

    const ctx = document.getElementById("tankCapacityChart");
    if (!ctx) return;

    if (charts.tankCap) charts.tankCap.destroy();

    charts.tankCap = new Chart(ctx, {
      type: "bar",
      data: {
        labels: ["Cylinder", "Total", "Oxidizer Used"],
        datasets: [{
          label: "Volume (m³)",
          data: [V_cyl, V_total, oxVol],
          backgroundColor: [
            "rgba(100, 217, 255, 0.8)",
            "rgba(91, 136, 255, 0.8)",
            "rgba(100, 217, 255, 0.5)"
          ],
          borderColor: "#64d9ff",
          borderWidth: 2,
          borderRadius: 8
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: true,
        indexAxis: "y",
        plugins: {
          legend: { display: false },
          tooltip: {
            backgroundColor: "rgba(15, 20, 40, 0.9)",
            titleColor: "#64d9ff",
            bodyColor: "#e8eef7",
            borderColor: "#64d9ff",
            borderWidth: 1
          }
        },
        scales: {
          x: {
            ticks: { color: "#b0b8d8" },
            grid: { color: "rgba(100, 217, 255, 0.1)" }
          },
          y: {
            ticks: { color: "#b0b8d8" }
          }
        }
      }
    });
  }

  // Update visualizations when compute buttons are clicked
  document.getElementById("engineComputeBtn")?.addEventListener("click", function() {
    setTimeout(updateVisualizations, 100);
  });

  document.getElementById("injComputeBtn")?.addEventListener("click", function() {
    setTimeout(updateVisualizations, 100);
  });

  document.getElementById("nozComputeBtn")?.addEventListener("click", function() {
    setTimeout(updateVisualizations, 100);
  });

  document.getElementById("fuelComputeBtn")?.addEventListener("click", function() {
    setTimeout(updateVisualizations, 100);
  });

  document.getElementById("grainComputeBtn")?.addEventListener("click", function() {
    setTimeout(updateVisualizations, 100);
  });

  document.getElementById("tankComputeBtn")?.addEventListener("click", function() {
    setTimeout(updateVisualizations, 100);
  });

  // Trigger reset all on page load to initialize the state
  resetAll();
});
