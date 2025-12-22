# Structure of FluxPrediction Class

This class represents a flux prediction results, e.g., from Flux Balance
Analysis or derived methods.

## Slots

- `algorithm`:

  Algorithm used for flux prediction.

- `ok`:

  The LP solver's code for the type of optimization result.

- `ok_term`:

  The LP solver's term for the type of optimization result.

- `stat`:

  Generic status (integer code) of the current basic solution of the
  optimization problem.

- `stat_term`:

  Generic status (character term) of the current basic solution of the
  optimization problem.

- `obj`:

  Objective value.

- `obj_sec`:

  Value of secondary objective function. E.g.: Summed absolute fluxes in
  [pfba](https://waschina.github.io/cobrar/reference/pfba.md).

- `fluxes`:

  Predicted flux values.

- `redCosts`:

  Predicted reduced costs (or "*dual value*") for reactions.
