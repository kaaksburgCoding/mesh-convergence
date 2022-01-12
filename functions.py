# -*- coding: utf-8 -*-

# Copyright 2013 SuperGen UK Centre for Marine Energy Research
# Copyright 2017-2020 Mathew Topper
#
# This file is part of convergence.
#
#    convergence is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    convergence is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with convergence.  If not, see <https://www.gnu.org/licenses/>.

"""
 -- convergence package  -------------------------------------------------------
 
   Performs several verification calculations given a file of grid spacings 
   and some observed quantity corresponding to each grid spacing.
 
   Computes:
   - order of convergence
   - Richardson extrapolation to zero grid spacing
   - grid convergence indices (GCI)
 
 --------------------------------------------------------------------------
 
Adapted from:
    
    NPARC Alliance CFD Verification and Validation Web Site
    Examining Spatial (Grid) Convergence
    verify.f90

    URL: http://www.grc.nasa.gov/WWW/wind/valid/tutorial/spatconv.html
    
    Nov '11: Updated to reflect Celik et al 2008.
"""

# Import built-in modules
from math import log


def order_of_convergence (value_1, value_2, value_3, ratio_21, ratio_32, 
                          omega=0.5, tol=1.E-4, max_iter=1e6): #, start_p=1.):
    
    """ Calculate the order of convergence values generated with three
    grids of reducing resolution (ie grid_1 is finest). The values of the grids
    are needed along with the ratios between them.
    
    An iterative method with under-relaxation is used to calculate the order
    of convergence as the refinement ratio is not necessarily constant.
    
    This has been modified to the method of Celik (2008).
    """
    
    # Set a maximum residual
    max_res = 1.E6

    print(value_1)
    print(value_2)
    print(value_3)

    # calculate the epsilons.
    epsilon32 = float(value_3 - value_2)
    epsilon21 = float(value_2 - value_1)

    print(epsilon21)
    print(epsilon32)
    
    # Calculate the fraction
    epfrac = epsilon32/epsilon21

    print('epfrac')
    print(epfrac)
    
    # Get the signed unit, s
    s = epfrac / abs(epfrac)
    
    # Initial guess at order of convergence, p
    p1 = (1. / log(ratio_21)) * abs(log(abs(epfrac))) # start_p
#p1 = (1. / log(ratio_21)) * abs(log(abs(epfrac)) + q)

    print("p0")
    print(p1)
    
    # Initialise the residual and number of iterations
    residual = 1.
    iterations = 0

    #as long as improvement of iterations (residual) is greater tol=1e-4
    while abs(residual) > tol:
        
        # Break if it's all gone bad
        if iterations > max_iter or residual > max_res:
            raise RuntimeError('Residual out of range or too many iterations')
        
        #p0 is initialised with the value of p1 to insert it into formula for q
        p0 = p1

        # Calculate q
        q = log((ratio_21 ** p0 - s) / (ratio_32 ** p0 - s))

        print("q")
        print(q)

        print(abs(log(abs(epfrac)) + q))

        # Calculate the pnew (whole formula including q)
        pnew = (1. / log(ratio_21)) * abs(log(abs(epfrac)) + q)
      #  p1 = (1. / log(ratio_21)) * abs(log(abs(epfrac)) + q)


        print("p1")
        print(p1)
        
        # Smooth the iteration relaxation step.
        # omega=0.5: p1=0.5p0+0.5pnew
        p1 = (1. - omega) * p0 + omega * pnew

        #caluclate difference (residual) between p of former iteration and
        #p of latest iteration
        residual = p1 - p0
        
        iterations += 1
    
    return p1


def richardson_extrapolate(value_j, value_k, ratio_jk, p):
    
    """ Estimate the zero grid spacing value using richardsons extrapolation and
    two grids of reducing resolution (ie grid_1 is finest). The refinement ratio
    is needed. The order of convergence, p, is also required.
    """
    
    f_exact = ( ratio_jk**p * value_j - value_k ) / ( ratio_jk**p - 1.0 )
    
    return f_exact


def error_estimates(value_1, value_2, f_exact):
    
    """ This routine returns the relative error and extrapolated 
    relative error. The values of the grids and needed along 
    with the extrapolated value.
    """
    
    # Get the approximate relative error
    e21a =  abs( (value_1 - value_2) / value_1 )
    
    # Get the extrapolated relative error
    e21ext = abs( ( f_exact - value_1 ) / f_exact )
    
    return e21a, e21ext


def gci(ratio_21, e21_approx, p):
    
    """ Calculate the fine and coarse grid convergence index for two grids of 
    reducing resolution (ie grid_1 is finest). The refinement ration between the 
    grids is required along with the approximate relative error (e21_approx) and 
    the order of convergence, p.
    """
    
    # Using a fixed safety factor as per Celik (2008)
    safety_factor = 1.25
    
    # Calculate the gci
    gci_fine = safety_factor * e21_approx / (ratio_21**p - 1.0)
    
    gci_coarse = ratio_21**p * gci_fine
    
    return gci_fine, gci_coarse


def asymptotic_ratio(gci_fine_21, gci_fine_32, ratio_21, p):
    
    """ Calculate the ratio in succesive Eps as defined at the bottom of page
    129 in Roache. If the ration is close to one then the asymptotic range has
    been reached.
    """
    
    ratio = ratio_21**p * ( gci_fine_21 / gci_fine_32)
    
    return ratio
