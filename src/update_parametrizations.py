import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description='Parse parametrization from SMASH analysis output')
parser.add_argument('source', help='File with SMASH analysis cross section output')
parser.add_argument('name', help='Prefix for C++ variable names')
args = parser.parse_args()

source = args.source
name = args.name

with open(source, 'r') as f:
    reaction = f.readline().lstrip('#initial ').rstrip('\n')
    for _ in range(3):
        f.readline()
    d = np.loadtxt(f)

sqrts = d[:,0]
elastic_contribution = d[:,2]

# Close to the production threshold, there can be numerical issues when doing
# the interpolation. To avoid these, it suffices to slightly modify the first
# energy.
sqrts[0] += 0.0025

def join_values(values, max_values_per_line, padding, precision):
    # center values on dot
    template = '{{left:>{0}}}.{{right:<{0}}}'.format(padding)
    formatted_values = []
    for value in values:
        value = "{{:.{}f}}".format(precision).format(value)
        left, right = value.split('.') if '.' in value else (value, '0')
        formatted_value = template.format(left=left, right=right)
        formatted_values.append(formatted_value)

    # arange values in table
    indent = '  '
    values_per_line = 0
    l = [indent]
    for value in formatted_values:
        l.append(value)
        if values_per_line < max_values_per_line:
            l.append(', ')
            values_per_line += 1
        else:
            l.append(',\n')
            l.append(indent)
            values_per_line = 0
    # avoid whitespaces at end
    if l[-1] == ', ':
        l[-1] = ''
    return ''.join(l)

s = '''/// Center-of-mass energy.
const std::initializer_list<double> {name}_RES_SQRTS = {{
{sqrts}
}};
/// Elastic {reaction} cross section contributions from decays.
///
/// These need to be subtracted from the interpolation of the PDG data on
/// elastic cross sections. This data was generated using the SMASH analysis
/// suite and should be updated when strange resonances are changed or added.
const std::initializer_list<double>{name} _RES_SIG = {{{elastic_contribution}}};
'''.format(name=name, reaction=reaction,
           sqrts=join_values(sqrts, 9, 0, 4),
           elastic_contribution=join_values(elastic_contribution, 4, 2, 11))

print s.rstrip('\n')
