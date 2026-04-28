# predict.brma Newdata Outcome Requirements

`predict.brma()` accepts `newdata` containing only the variables used by the
requested prediction type.

## Required Outcome Columns

- `type = "terms"` / `type = "terms.scale"`: no outcome columns are required,
  except normal PET/PEESE predictions with `bias_adjusted = FALSE`, which need
  `sei` or `vi` for the bias term.
- `type = "estimate"` for new data: no observed outcome is required; normal
  PET/PEESE predictions with `bias_adjusted = FALSE` still need `sei` or `vi`.
- `type = "response"`: sampling-distribution inputs are required.
  - Normal models need `sei` or `vi`; `yi` is not required.
  - Binomial GLMMs need `n1i` and `n2i`; observed `ai`/`ci` are not required.
  - Poisson GLMMs need `t1i` and `t2i`; observed `x1i`/`x2i` are not required.

Unused outcome columns are filled internally with dummy values so the shared
data parser can still construct a complete `RoBMA_data` object.
