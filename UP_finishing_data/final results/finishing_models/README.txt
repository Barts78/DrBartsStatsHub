UNDER PRESSURE: LONGITUDINAL FINISHING MODELS

Input: /Users/Barts78/Documents/Dropbox files/AI lfc data/R files/DrBartsHubPage/DrBartsStatsHub/UP_finishing_data/final results/finishing_five_seasons/player_seasons.csv
Min shots per qualifying season: 30
Min qualifying seasons per player: 2
Main cohort: 152 players; 444 player-seasons

DEFINITION: FVA per shot = sum(xGOT - xG) / number of usable non-penalty shots.
This script uses the FVA values already computed by the earlier pipeline; it does not reconstruct missing shot-level values.
Unknown player ID -1 excluded to prevent merging unrelated players.
Within-player slope: weighted fixed-effects trend across calendar seasons.
Uncertainty: player-cluster bootstrap percentile intervals (resample players, preserving each player trajectory).
Shot volume is a precision proxy; heterogeneity in shot context and model calibration remains unadjusted.
Year-to-year correlations pool adjacent qualifying seasons and are descriptive, not reliability-corrected.
Complete-xGOT sensitivity excludes any player-season with missing xGOT; may select non-representative data.
The goals threshold is intentionally not used to define the main cohort (selection on outcome).
No age or position inference is made unless a trustworthy age column is supplied.
Calendar-year trends can reflect aging, tactical changes, team changes, selection and measurement noise.
Confidence intervals excluding zero are not evidence that coaching caused change.

OUTPUTS: within_player_trends.csv; adjacent_season_pairs.csv; year_to_year_correlations.csv;
player_trends_descriptive.csv; qualifying_player_seasons.csv; cohort_season_summary.csv; figures.

NEXT: add position, date of birth / age, and shot-level covariates for a fuller hierarchical model.
