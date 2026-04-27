**brma class**
cluster/estimate follow-up
  - GLMM cluster likelihood: validate against high-node/brute-force quadrature.
  - selection/weighted cluster likelihood: document conditional likelihood semantics; no full cluster-selection process yet.
  - cluster_likelihood.n_gamma: stress-test default and numerical robustness.

missing features
  - vif 	Variance inflation factors (search how GIF works in metafor)
  - vcov	Variance-covariance matrix of coefficients (requires more methodological development)
  - robust	Cluster-robust inference (requires more methodological development)
  - baujat	Baujat plot (skip too computationally difficutl)
  - gosh	GOSH analysis

missing wrappers:
  - confint	credible intervals (for all model parameters i.e., merges output from summary tables )
  - fitted	Fitted values (currently only via predict)
  - ranef	Extract random effects (should this be the offsets? see metafor)

extensions of existing methods
  - add cluster argument to loo / residuals / forest plot

plot_weightfunction (should be in BayesTools)
Expand the function to use `show_data` argument. if TRUE, the function should show p-values as jitter on the x-axis
- compute p-values by transformation of yi/sei
- note that re-scale options spaces the cutpoints equally (not mapping to p-values directly)
- note that two-sided only weightfunction with bselmodel use two sided p-values, other use one-sided p-values

prompts:
Now I want you to fully implement the plan. Take your time, use subagents. Work like a senior developer.
Remmember that the unit tests are currently broken so you cannot relly on them. Once you are done, we will go over them and fully test the remaining steps.
