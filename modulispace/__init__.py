"""
This package implements moduli spaces over a K3 surface X with Pic X = Z.

The functionality of the package comprises:
    - Computation of the Movable cone.
    - Computation of Flopping contractions (spherical and sum flop)
    - Printing the computed wall-and-chamber decomposition of the movable cone.
    - Showing plots of the walls in the alpha-beta-plane.

The functions read_walls_file, write_walls_file have the functionality to
write/read the computed walls to/from a file.

Algorithm: This package implements the wall-and-chamber decomposition from
Bayer--Macrí: "MMP for moduli of sheaves on K3s via wall-crossing"
in Inventiones mathematicae 198 (3), 505-590.

Overview of algorithm from loc.cit.:

X K3 surface.

Pic X = Z H.

H_alg = Z + ZH + Z with (a,k,b).(c,k',d) = kk'-ad-bc.

M = M(r,D,c) the moduli space with Mukai vector v = (r,D,c).

Then Pic M = (r,D,c)^perp.

Mov M is cut out in Pos M by:
    a^ perp cap v^perp
for all a with
    - a^2 = -2 or  (Brill-Noether contraction)
    - a^2 = 0 and v.a = 1 or 2.  (Hilbert-Chow or Li-Gieseker-Uhlenbeck contraction)

The flopping contractions are given by
    a^perp cap v^perp
for all a such that
    - a^2 = -2 and 0 < v.a <= v^2/2, or
    - v = a + b for a,b two positive classes in the sense of loc.cit.
"""


from .modulispace import ModuliSpace
from .modulispace import read_walls_file, write_walls_file
from .mukaivector import mVec
from .k3surface import K3Surface
from .walls import Wall



__all__ = ["ModuliSpace", "mVec", "K3Surface", "Wall", "read_walls_file", "write_walls_file", ]