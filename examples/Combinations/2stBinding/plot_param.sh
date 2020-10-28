#!/bin/sh

# Display chemical shift difference between the two states
chemex plot_param -p Output/STEP2/All/Parameters/fitted.toml -n DW_AB

# Display R2 rates for both states
chemex plot_param -p Output/STEP2/All/Parameters/fitted.toml -n R2
