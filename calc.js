// calc.js -- ESKD risk from the exported Fine-Gray models (model.js). Shared by index.html and test_calc.js.
(function (root) {
  "use strict";

  // natural-spline basis from exact piecewise polynomials (see export_model.R)
  function splineBasis(sp, x) {
    var br = sp.breaks, last = br.length - 1, rows, u;
    if (x < br[0]) { rows = sp.left; u = x - br[0]; }
    else if (x > br[last]) { rows = sp.right; u = x - br[last]; }
    else {
      var i = 0; while (i < last - 1 && x >= br[i + 1]) i++;
      rows = sp.inner[i]; u = x - br[i];
    }
    return rows.map(function (c) { var v = 0, p = 1; for (var k = 0; k < c.length; k++) { v += c[k] * p; p *= u; } return v; });
  }

  // design row in the model's column order: linear terms, then each spline's basis columns
  function design(m, x) {
    var row = m.lin.map(function (v) { return x[v]; });
    m.spl.forEach(function (v) { row = row.concat(splineBasis(m.splines[v], x[v])); });
    return row;
  }
  function linpred(m, x) {
    var row = design(m, x), lp = 0;
    for (var j = 0; j < row.length; j++) lp += (row[j] - (m.xbar ? m.xbar[j] : 0)) * m.coef[j];
    return lp;
  }
  function gridIndex(M, t) { return Math.min(Math.round(t / M.grid_step), M.eskd.H0.length - 1); }

  // cumulative incidence at time t (years, on the 0.05 grid); same form for Fine-Gray and the death-censored Cox
  function risk(M, which, x, t) {
    var m = M[which], h = m.H0[gridIndex(M, t)];
    return 1 - Math.exp(-h * Math.exp(linpred(m, x)));
  }
  function curve(M, which, x) {
    var e = Math.exp(linpred(M[which], x));
    return M[which].H0.map(function (h) { return 1 - Math.exp(-h * e); });
  }
  // x for the spline models from clinical units
  function inputs(age, egfr, utpMg, albumin) {
    return { Age: age, eGFR: egfr, log_utp: utpMg == null ? null : Math.log1p(utpMg), Albumin: albumin };
  }
  function palladiniStage(egfr, utpMg) { return ["I", "II", "III"][(egfr < 50 ? 1 : 0) + (utpMg >= 5000 ? 1 : 0)]; }

  var api = { splineBasis: splineBasis, risk: risk, curve: curve, inputs: inputs, palladiniStage: palladiniStage };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.ESKD = api;
})(this);
