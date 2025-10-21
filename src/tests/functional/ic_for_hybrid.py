#!/usr/bin/env python3

# ===================================================
#
#    Copyright (c) 2024
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
# ===================================================

import argparse
import sys
import subprocess as subp
import numpy as np
import pandas as pd


def cartesian_to_milne(t, z):
    tau = np.sqrt(t**2 - z**2)
    eta = 0.5*np.log((t+z)/(t-z))
    return tau, eta


parser = argparse.ArgumentParser(
    description='Functional test for iso-tau initial conditions.')
parser.add_argument('--source', required=True, help='SMASH top level')
parser.add_argument('--binary', required=True, help='build directory')
parser.add_argument('--skip-smash', action='store_true')
args = parser.parse_args()

# Each functional test should have its own output folder as functional tests
# might be run in parallel and, in that case, smash would fail if different
# smash executables would try to print to the same output folder.
output_directory = args.binary + "/functional_test_output/ic_for_hybrid"
# Run smash with appropriate configuration
if not args.skip_smash:
    smash_executable = args.binary + "/smash"
    smash_config_file = args.source + "/input/config.yaml"
    smash_config_output = \
        "Output: {Initial_Conditions: {Format: [for_vHLLE, Binary, Oscar2013]}}"
    smash_config_ic = \
        "Modi: {Collider: {Initial_Conditions: {Type: Constant_Tau}}}"
    run = subp.run([smash_executable, "-i", smash_config_file,
                    "-c", smash_config_output,
                    "-c", smash_config_ic,
                    "-o", output_directory, "-f", "-q"],
                   cwd=args.binary)
    try:
        run.check_returncode()
    except subp.CalledProcessError:
        print("SMASH did not run properly. Please check that " +
              output_directory + " is clean.")
        sys.exit(1)

# Check if output is consistent
ascii_file = output_directory + '/SMASH_IC_for_vHLLE.dat'
oscar_file = output_directory + '/SMASH_IC.oscar'

# Read output files
with open(ascii_file, 'r') as f:
    next(f)  # ignore first line
    ascii_quantities = f.readline().split()[1:]
with open(oscar_file, 'r') as f:
    oscar_quantities = f.readline().split()[2:]
ascii = pd.read_csv(ascii_file, comment='#', delimiter=' ',
                    names=ascii_quantities, index_col=False)
oscar = pd.read_csv(oscar_file, comment='#', delimiter=' ',
                    names=oscar_quantities, index_col=False)

# Oscar particles are the ascii particles plus spectators
number_of_spectators = len(oscar) - len(ascii)
if (number_of_spectators < 0):
    print("output for vHLLE should not have more particles than OSCAR.")
    sys.exit(1)
oscar['tau'], oscar['eta'] = cartesian_to_milne(oscar['t'], oscar['z'])
oscar['mt'], oscar['Rap'] = cartesian_to_milne(oscar['p0'], oscar['pz'])

quantities_to_compare_exact = ['x', 'y', 'pdg']
# px and py are given with default precision in "for_vHLLE", but with %.9 in OSCAR.
# If this changes, they should be compared exactly
quantities_to_compare_fuzzy = {
    'px': 1.e-5, 'py': 1.e-5, 'eta': 1.e-4, 'tau': 1.e-5, 'mt': 1.e-5, 'Rap': 1.e-5}
comparison_exact = pd.merge(ascii[quantities_to_compare_exact],
                            oscar[quantities_to_compare_exact],
                            on=quantities_to_compare_exact, how='left', indicator='exists')
if set(comparison_exact['exists']) != {'both'}:
    print("Failed comparison of quantities that should be exactly equal.")
    sys.exit(1)

# - Remove spectators from oscar and check if they are nucleons. Since there
# is no ID in the ascii output, this has to be done by comparing the quantities
# we know to be exact.
spectators = pd.concat([ascii[quantities_to_compare_exact],
                        oscar[quantities_to_compare_exact]]).drop_duplicates(keep=False)
if len(spectators) != number_of_spectators:
    sys.exit(1)
if not set(spectators['pdg']).issubset([2212, 2112]):
    print("Some spectators seem to not be nucleons.")
    sys.exit(1)
oscar = oscar.drop(spectators.index.values).reset_index()

# The required precision is arbitrary.
for quantity, tolerance in quantities_to_compare_fuzzy.items():
    comparison_fuzzy = abs((ascii[quantity] - oscar[quantity])/ascii[quantity])
    if any(comparison_fuzzy > tolerance):
        print("Relative difference in " + quantity + " above " + tolerance)
        sys.exit(1)
