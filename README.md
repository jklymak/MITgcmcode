Energy Diagnostics
------------------

WARNING: There appear to be some issues with `FIND_RHO` in here that need to be worked on for modern MITgcm versions.

Need:

* `energy_diagnostics_fill.F`
* `do_statevars_diags.F`

in your `code` directory

You also need to add a `data.diagnostics` file like the example here. And `packages.conf` should have `diagnostics` turned on.

## Dissipation due to overturns 

WARNING: for "newer" MITgcm versions, just use `pkg/kl10`

This is based off Klymak and Legg 2010, Ocean Modelling.

You need

* `PP81.h`
* `PP81_DIAGS.h`
* `pp81_calc.F`
* `pp81_output.F`
* `pp81_readparms.F`

in your `code` directory.  In order to get the dissipation output, you need to edit `data.pp81` in your `input` directory.



