"""
walls.py implements walls in the sense of Bayer--Macrí.

It adds the functionality of writing and reading walls from a file.


The file format is as follows (with # annotation and explanation):

4   # Degree of K3

(1, 0, -9)  # Mukai vector of Moduli space

# Then each line is a new wall as follows: For example

(-3, 7, -33);4;(7, -15, 63);None;7

Seperated by a semicolon we have:

# Mukai vector inducing the wall;

# Type of the wall (see below in the class Wall);

# induced movable divisor;

# contracted divisor (if not flopping contraction), else None

# codim (if flopping contraction), else None

"""
from __future__ import annotations
from modulispace.mukaivector import mVec
from typing import Optional


class Wall:
    """
    Class corresponding to a wall in a Modulispace over a k3 surface.
    Type of the wall can be any of [NON, BN, HC, LGU, LAGRANGIAN, FLOP_SPHERE, FL_SUM].
    """
    NON = -1  # Not a Wall
    BN = 0  # Brill--Noether
    HC = 1  # Hilbert-Chow
    LGU = 2  # Li-Gieseker-Uhlenbeck
    LAGRANGIAN = 3
    FLOP_SPHERE = 4  # Flop sum
    FL_SUM = 5  # Flop sperical
    dict = {NON: "No Contraction", BN: "Brill-Noether", HC: "Hilbert-Chow",
            LGU: "Li-Gieseker-Uhlenbeck", FLOP_SPHERE: "Spherical Flop",
            FL_SUM: "Sum Flop", LAGRANGIAN: "Lagrangian Fibration"}
    dict_short = {NON: "No Contraction", BN: "BN", HC: "HC",
                  LGU: "LGU", FLOP_SPHERE: "Flop",
                  FL_SUM: "Flop", LAGRANGIAN: "Lagrangian"}

    def __init__(self, special_class: mVec, wall_type: int, induced_movable_divisor: mVec,
                 contracted_divisor: Optional[mVec] = None, codim: int = None):
        """

        :param special_class: class inducing the wall
        :param wall_type: type of the wall [NON, BN, HC, LGU, LAGRANGIAN, FLOP_SPHERE, FL_SUM].
        :param induced_movable_divisor: the movable divisor inducing the contraction
        :param contracted_divisor: the contracted divisor (if it exists)
        :param codim: codimension of the contracted locus
        """

        self.special_class = special_class
        self.wall_type = wall_type
        self.induced_divisor = induced_movable_divisor
        self.contracted_divisor = contracted_divisor
        self.codim = codim
        if (wall_type == -1 or wall_type >= 4) and contracted_divisor is not None:
            raise ValueError

    def get_vec(self) -> mVec:
        return self.special_class

    def get_movable_divisor(self) -> mVec:
        return self.induced_divisor

    def get_contracted_divisor(self) -> mVec:
        return self.contracted_divisor

    def get_type(self) -> int:
        return self.wall_type

    def set_type(self, wall_type: int) -> None:
        self.wall_type = wall_type

    def get_codim(self) -> int:
        return self.codim

    def get_type_str(self) -> str:
        return self.dict[self.wall_type]

    def get_type_short_str(self) -> str:
        return self.dict_short[self.wall_type]

    def make_primitive(self) -> Wall:
        """
        Changes the induced movable divisor/contracted divisor to its primitive  integral divisor.
        """
        self.induced_divisor = self.induced_divisor.get_primitive()
        if self.contracted_divisor is not None:
            self.contracted_divisor = self.contracted_divisor.get_primitive()
        return self

    def __eq__(self, other: Wall) -> bool:
        return self.induced_divisor == other.induced_divisor

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        string = self.dict[self.wall_type]
        if self.wall_type == self.FL_SUM or self.wall_type == self.FLOP_SPHERE:
            string += " of codimension " + str(self.codim)
        string += "\n"
        string += "Vector inducing Wall: " + str(self.special_class) + "\n"
        string += "Movable Divisor:      " + str(self.induced_divisor) + "\n"
        if self.contracted_divisor is not None:
            string += "Contracted Divisor:   " + str(self.contracted_divisor) + "\n"
        return string
