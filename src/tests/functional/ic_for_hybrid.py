#!/usr/bin/env python3
#===================================================
#
#    Copyright (c) 2024
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
#===================================================

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

output_dir = "./test_output/ic_for_hybrid"
# Run smash with apropriate configuration
if not args.skip_smash:
    smash_config_file = args.source+"/input/config.yaml"
    smash_config_output = \
        "Output: {Initial_Conditions: {Format: [ASCII, Binary, Oscar2013]}}"
    smash_config_ic = \
        "Modi: {Collider: {Initial_Conditions: {Type: Constant_Tau}}}"
    run = subp.run(["smash", "-i", smash_config_file,
            "-c", smash_config_output,
            "-c", smash_config_ic,
            "-o", output_dir,"-f","-q"],
            cwd=args.binary)
    try:
        run.check_returncode()
    except subp.CalledProcessError:
        print("SMASH did not run properly. Please check that " + output_dir +
              " is clean.")
        sys.exit(1)

# Check if output is consistent
ascii_file = args.binary + '/' + output_dir + '/SMASH_IC.dat'
oscar_file = args.binary + '/' + output_dir + '/SMASH_IC.oscar'

# Read output files
with open(ascii_file, 'r') as f:
    next(f) #ignore first line
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
    print("ASCII should not have more particles than OSCAR.")
    sys.exit(1)
oscar['tau'], oscar['eta'] = cartesian_to_milne(oscar['t'], oscar['z'])
oscar['mt'], oscar['Rap'] = cartesian_to_milne(oscar['p0'], oscar['pz'])

quantities_to_compare_exact = ['x','y','pdg']
# px and py are given with default precision in ASCII, but with %.9 in OSCAR.
# If this changes, they should be compared exactly
quantities_to_compare_fuzzy = ['px','py','eta','tau','mt','Rap']
compare_df = pd.merge(ascii[quantities_to_compare_exact],
                  oscar[quantities_to_compare_exact],
                  on=quantities_to_compare_exact, how='left', indicator='exists')
if set(compare_df['exists']) != {'both'}:
    print("Failed comparison of quantities that should be exactly equal.")
    sys.exit(1)

# - Remove spectators from oscar and check if they are nucleons. Since there
# is no ID in the ascii output, this has to be done by comparing the quantities
# we know to be exact.
spectators = pd.concat([ascii[quantities_to_compare_exact],
                        oscar[quantities_to_compare_exact]]).drop_duplicates(keep=False)
if len(spectators) != number_of_spectators:
    sys.exit(1)
if not set(spectators['pdg']).issubset([2212,2112]):
    print("Some spectators seem to not be nucleons.")
    sys.exit(1)
oscar = oscar.drop(spectators.index.values).reset_index()

# The required precision is arbitrary.
for quantity in quantities_to_compare_fuzzy:
    if any(abs((ascii[quantity] - oscar[quantity])/ascii[quantity]) > 10**-5):
        print("Relative difference in "+quantity+" above 10⁻⁵.")
        sys.exit(1)
