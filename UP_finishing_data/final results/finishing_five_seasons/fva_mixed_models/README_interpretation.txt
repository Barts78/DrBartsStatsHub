Under Pressure | Mixed-effects finishing analysis
Primary model: player-season FVA/shot ~ season fixed effects + player random intercept.
Comparator: xG/shot with same design and shot-count precision proxy.
Weights: usable shots / median usable shots for FVA; shots / median shots for xG.
ICC = between-player variance / (between-player variance + residual variance at median shot count).
Weighted models have shot-count-dependent residual variance; ICC is NOT a universal property.
Bootstrap: parametric, conditional on model assumptions. Small ICC estimates may be boundary/singular.
No ages or positions in source data; cannot infer coaching effects or age curves.
Missing xGOT sensitivity changes the cohort; do not interpret as causal effect of missingness.
A large xG/shot ICC can reflect stable positions, tactics and chance types.
Player random intercepts are shrinkage estimates, not rankings of true finishing talent.
Selection: only player-seasons above threshold; survivor and role bias remain.
Model is an exploratory Gaussian approximation to aggregated per-shot FVA; a shot-level hierarchical model would be stronger.
