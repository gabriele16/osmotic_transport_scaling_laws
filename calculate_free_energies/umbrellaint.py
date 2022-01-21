#!/usr/bin/env python

import re
import warnings
import argparse
from scipy.constants import Boltzmann as KB
from scipy.constants import Avogadro as NA
import numpy as np
#import pandas as pd
from argparse import RawTextHelpFormatter
from scipy.integrate import trapz
import sys

description = """Umbrella Integration.
The implementation follows "Kaestner, J. (2011). Umbrella sampling. 
Wiley Interdisciplinary Reviews: Computational Molecular Science, 1(6), 932-942."
### metafile format:
/window/data window_center sprint_konst [Temperature]
### window data file format:
time_step coordinate (1-dimentional)
"""

arg_parser = argparse.ArgumentParser(
    description=description, formatter_class=RawTextHelpFormatter)
arg_parser.add_argument('-o', '--output',
                        metavar='Output free energy file',
                        default='free_energy.dat', dest='out_put',
                        help="Optional, use 'free_py.txt' as default",)
arg_parser.add_argument('-T', '--temperature',
                        metavar='Temperature',
                        dest='temperature', default=-1, type=float,
                        help="Optional, set a default temperature.")
arg_parser.add_argument('-c', '--confidence',
                        default=2, type=int, metavar='1|2|3', dest='confidence',
                        choices=[1, 2, 3 ],
                        help='confidence interval for error calculation in terms of sigma, i.e. 1*sigma, 2*sigma or 3*sigma')
arg_parser.add_argument('-X', '--range',
                        nargs=2, default=None, metavar='Range of xi',
                        dest='range', type=float,
                        help="Range of reaction coordinate.")
arg_parser.add_argument('meta_file',
                        nargs=None,
                        help='Meta file name')
arg_parser.add_argument('max_bin',
                        type=int, nargs=None,
                        help='How many bins were used in integration.')
arg_parser.add_argument('start',
                        type=int, nargs=None, default = 5000,
                        help='get data after start lines')
arg_parser.add_argument('end',
                        type=int, nargs=None, default = 55000,
                        help='finish analysis at end lines')
arg_parser.add_argument('-s','--split',
                        type=str, nargs=None, dest ='split',default='no',
                        choices=['no','upper_half','lower_half'], metavar='no|upper_half|lower_half',
                        help='split data in two or not')



args = arg_parser.parse_args()
alvars = vars(args)

#convert Bohr to nanometer
#use in xi, as the CV is printed in Bohr as default in CP2K
bohr2nm=0.0529177

# Variables

temperature = alvars['temperature']
max_bin = alvars['max_bin']  # how many bins were used in integration
xi_range = alvars['range']
out_put_file = open(alvars['out_put'], 'w')
out_put_file.write('#r PMF MF\n')
out_put_file.close()
out_put_file = open(alvars['out_put'], 'a')
start = alvars['start']
end = alvars['end']
confidence = float(alvars['confidence'])
split=alvars['split']

kb = 0.0083144621
if xi_range:
    print(xi_range)
    assert xi_range[0] < xi_range[1], "Give the right range!"

meta_file = open(alvars['meta_file'], 'r')
window_info, minimum = [], []


class NoTemperatureError(Exception):
    r"""No temperature error."""

    pass


for line in meta_file:
    if not re.search('^#', line) is None:
        continue
    line = re.split('\s+', line.strip())
    if temperature == -1 and len(_line) != 4:
        raise NoTemperatureError("You have not set temperature for this "
                                 "window or a global temperature!")
    i = 0
    cv_vs_t=[]
    for l in open(line[0]):
       i+=1
       columns = l.split()
       cv_instantaneous = bohr2nm * float(columns[1])
       if i >= start and i < end:
          cv_vs_t.append(cv_instantaneous)
       elif i >= end:
          break
    
    window_data = np.array(cv_vs_t)

    if split == "upper_half":

        window_data = window_data[int(len(window_data) /2) : ]

    elif split == "lower_half":  

        window_data = window_data[: int(len(window_data) /2)  ]

        
    center = float(line[1])
    spring_konst = float(line[2])
    kbT = kb * temperature
    window_info.append([window_data.mean(), window_data.var(),
                         center, spring_konst, kbT, len(window_data)])
    minimum.append(window_data.min())
    minimum.append(window_data.max())

window_info = np.array(window_info)
window_info = window_info[np.argsort(window_info.T[2])]

if xi_range:
    if min(minimum) < xi_range[0] or max(minimum) > xi_range[1]:
        warnings.warn("Warning, xi range exceeds the sample range!",
                      UserWarning)

xi_range = xi_range or [min(minimum), max(minimum)]
xis = np.linspace(xi_range[0], xi_range[1], max_bin)
xi_mean_w = window_info.T[0][:, np.newaxis]
xi_var_w = window_info.T[1][:, np.newaxis]
xi_center_w = window_info.T[2][:, np.newaxis]
k_w = window_info.T[3][:, np.newaxis]
kbT_w = window_info.T[4][:, np.newaxis]
Ni =window_info.T[5][:,np.newaxis]

pbc_xis = xis

# \partial A/\partial \xi_{bin} =
# \sum_i^{window} P_i(\xi_{bin})/(\sum_i^{window} P_i(\xi_{bin})) \times
# \partial A_i^u/\partial \xi_{bin}
# \partial A_i^u / \xi_{bin}, with shape (n_window, n_xi)
dAu_dxis = kbT_w * (pbc_xis - xi_mean_w) / xi_var_w -\
        k_w * (pbc_xis - xi_center_w)



var_dA_dxi = kbT_w**2 * ( 2.0* (pbc_xis - xi_mean_w)**2  + xi_var_w )
var_dA_dxi /= (Ni * (xi_var_w**2) )


# N_iP_i(\xi_{bin}), with shape (n_window, n_xi),
# all Nis are same in this case
pb_i = 1/np.sqrt(2 * np.pi) * 1 / np.sqrt(xi_var_w) *\
        np.exp(-0.5 * (pbc_xis - xi_mean_w) ** 2 / xi_var_w)

p_xi = (pb_i*Ni)/np.sum(pb_i*Ni,axis=0)


normalization = np.sum(p_xi, axis=0)

if (normalization.any() != 1.0):
   print("error weights do not normalize to 1")
   print(normalization)
   sys.exit()

# error calculation
p_xi_sq = p_xi*p_xi
var_dA_dx = np.sum(var_dA_dxi*p_xi*p_xi,axis=0)
mean_sigma_xi = np.mean(np.sqrt(xi_var_w))
var_dA = np.mean(var_dA_dx)
var_dA *= ( max(xis) - min(xis) *(  mean_sigma_xi * np.sqrt(2.0*np.pi) ) -2.0*mean_sigma_xi**2 )


# calculation of sum( dA/dxi *p_xi)  for all windows
dA_dxis = np.sum(dAu_dxis * p_xi, axis=0)

# error is taken as standard deviation * confidence interval
error_dA_vec = np.array([np.sqrt(var_dA) * confidence]*len(xis))

pmf = np.array([trapz(dA_dxis[xis <= r], xis[xis <= r]) for r in xis]) 
pmf -= np.min(pmf)

np.savetxt(out_put_file, np.vstack([xis, pmf, error_dA_vec ,dA_dxis]).T, fmt="%.6f")
out_put_file.close()
