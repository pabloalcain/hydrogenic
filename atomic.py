# -*- coding: utf-8 -*-
"""
The ideal set up of this is a class called 'Hydrogenic' that can be
instantiated for each type of atom. Unfortunately this cannot be done
easily because of the use of common blocks in the fortran back-end

"""

import numpy as np
import radial
global potential
potential = None

def set_potential(Z, Zs, alpha):

  """
  Set the potential for the simulation to Coulomb and Yukawa
  V(r) = Z/r + Zs/r*exp(-alpha*r)

  Parameters
  ----------
  Z : float
    Un-screened charge (in electron units) of the nucleus

  Zs: float
    Screened charge (in electron units) of the nucleus

  alpha: float
    Yukawa screening parameter, inverse of screening length

  Notes
  -----
    The charge is negative for protons in units of electron charge!
  """
  global potential
  radial.potential(Z, Zs, alpha)
  potential = lambda r: Z/r + Zs/r * np.exp(-alpha*r)


def bound(n, lk, eq, eps=1.0e-15):
  """
  Calculate the bound state of the Schrödinger or Dirac equation for
  the hydrogenic atom

  Parameters
  ----------
  n : integer
    Principal quantum number

  lk: integer
    Orbital angular momentum for Schrodinger, k coupling for Dirac

  eq: {'S', 'D'}
    Equation to use; 'S' stands for Schrödinger, 'D' for Dirac equation

  eps: float, optional
    Accuracy goal (the default is 1.0e-15)

  Returns
  -------
  E: float
    Energy of the bound state

  R: numpy array
    Positions of the interpolation grid

  P: numpy array
    Radial wavefunction u

  Q: numpy array
    Radial wavefunction derivative u'

  Notes
  -----

  `P` and `Q` *are not* the radial wavefunctions and derivative, but rather

  .. math::
    P(r) = u_{nl}(r) = R_{nl}(r)\,r\\
    Q(r) = u'_{nl}(r) = R'_{nl}(r)\,r + R_{nl}

  """
  if eq == 'S':
    E = radial.schroed_bound(n, lk, eps)
  elif eq == 'D':
    E = radial.dirac_bound(n, lk, eps)
  else:
    raise ValueError("eq should be 'D' [for Dirac] or 'S' [for Schrödinger]")
  r = radial.radwf.rad[:radial.radwf.ngp]
  p = radial.radwf.p[:radial.radwf.ngp]
  q = radial.radwf.q[:radial.radwf.ngp]
  return E, r, p, q

def free(e, lk, eq, eps=1.0e-15):
  """
  Calculate the free state of the Schrödinger or Dirac equation for the
  hydrogenic atom

  Parameters
  ----------
  e : float
    Energy of the free state

  lk: integer
    Orbital angular momentum for Schrodinger, k coupling for Dirac

  eq: {'S', 'D'}
    Equation to use; 'S' stands for Schrödinger, 'D' for Dirac equation

  eps: float, optional
    Accuracy goal (the default is 1.0e-15)

  Returns
  -------
  IPS: float
    Internal Phase Shift of the wave

  CPS: float
    Coulomb Phase Shift of the wave

  ETA: float
    ETA (?)

  R: numpy array
    Positions of the interpolation grid

  P: numpy array
    Radial wavefunction u

  Q: numpy array
    Radial wavefunction derivative u'

  Notes
  -----

  `P` and `Q` *are not* the radial wavefunctions and derivative, but rather

  .. math::
    P(r) = u_{nl}(r) = R_{nl}(r)\,r\\
    Q(r) = u'_{nl}(r) = R'_{nl}(r)\,r + R_{nl}

  """
  if eq == 'S':
    ips, cps, eta = radial.schroed_free(e, lk, eps)
  elif eq == 'D':
    ips, cps, eta = radial.dirac_free(e, lk, eps)
  else:
    raise ValueError("eq should be 'D' [for Dirac] or 'S' [for Schrödinger]")
  r = radial.radwf.rad[:radial.radwf.ngp]
  p = radial.radwf.p[:radial.radwf.ngp]
  q = radial.radwf.q[:radial.radwf.ngp]
  return ips, cps, eta, r, p, q
