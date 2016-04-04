#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Tiny example to show how to use the atomic module
"""

import atomic as at
import pylab as pl

at.set_potential(-1, 0, 0)
pl.figure()
for n, l in zip((1, 2, 2), (0, 0, 1)):
  E, r, p, q = at.bound(n, l, 'S')
  print "Energy of Schrödinger for n = {n} and {l} is {E}".format(n=n, l=l, E=E)
  pl.plot(r, p/r, label='n = {n}, l = {l}'.format(n=n, l=l))
pl.legend()
pl.xlabel('r')
pl.ylabel(r'$\Psi(r)$')
pl.title('Hydrogen bound states')
pl.savefig('h_bound.png')
pl.show()
