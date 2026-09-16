# Output semantics

`transition_ready_endpoint_support.csv` is the main hand-off file for the next transition stage. Each row is an actual HR parameter vector that passed either the v3.9 optimum or the near-optimal identifiability search. `group_support_weight` ensures equal cell contribution and prevents cells with more alternative solutions from receiving more biological weight.

`core_q75_secure` is a diagnostic flag only. A cell is not removed when this flag is false. It is true when rheobase, q=.75 firing rate and q=.75 mean ISI are each secured either by full parameter identifiability or by ≤20% variation across all evaluated near-optimal alternatives.

The q=.50 layer is extended/supporting. There is no imputation or extrapolation for the unsupported cell.

The support-member weights are not Bayesian probabilities and must not be interpreted as confidence probabilities.
