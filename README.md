# ESKD risk calculator for renal AL amyloidosis

Individual 1-, 3- and 5-year risk of end-stage kidney disease (dialysis or transplant) in patients with AL amyloidosis
and kidney involvement, accounting for the competing risk of death. Open it at
**https://shankara-a.github.io/eskd-calculator/**. It runs entirely in the browser; nothing entered leaves the page.

**Model.** Fine-Gray competing-risks model with age and natural splines (3 df) of eGFR and log(1 + 24-hour urine
protein, mg/day), fitted in 1,282 patients with renal AL amyloidosis evaluated at the Boston University Amyloidosis
Center, 1994-2023. Time zero is the first center visit, when the labs are drawn; death without ESKD is the competing
event. Internally validated (bootstrap and repeated 10-fold cross-validation over 20 multiply imputed datasets).
Not yet externally validated. A fallback model (age, eGFR, serum albumin) is used when urine protein is not available.

**Files.** `index.html` (page), `calc.js` (prediction), `model.js` (exported model: coefficients, spline bases as
piecewise polynomials, baseline cumulative hazards on a 0.05-year grid, Palladini stage summaries). Aggregate model
parameters only; no patient-level data. Predictions match the R implementation (`cmprsk::crr` / `predict.crr`) to
within 1e-14.

**Citation.** Anand S, Staron A, Burke N, Sanchorawala V, Verma A. An individualized competing-risks model for kidney
failure in renal AL amyloidosis. ASH Annual Meeting 2026, abstract abs26-10165.

For research use only; not a substitute for clinical judgement.
