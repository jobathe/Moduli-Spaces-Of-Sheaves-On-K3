"""
modulispace.py implements moduli spaces over a K3 surface X with Pic X = Z.

The functionality of the class Modulispace comprises:
    - Computation of the Movable cone.
    - Computation of Flopping contractions (spherical and sum flop)
    - Printing the computed wall-and-chamber decomposition of the movable cone.
    - Showing plots of the walls in the alpha-beta-plane.

The functions read_walls_file, write_walls_file have the functionality to
write/read the computed walls to/from a file.

Algorithm: This file implements the wall-and-chamber decomposition from
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
import math
import sys
from modulispace.mukaivector import mVec
from modulispace.qequation import QEquation
from modulispace.lattice import Lattice
from modulispace.walls import Wall
from modulispace.k3surface import K3Surface
from typing import Union, Optional, TextIO
try:
    from sympy import simplify, Matrix
    from sympy.core.power import integer_nthroot
    from sympy.ntheory.primetest import is_square
    from sympy.core.numbers import Rational
except ImportError:
    print("SymPy is required for this program. Please install it via \n pip install sympy\n")
    sys.exit()
try:
    import matplotlib.pyplot as plt
    disable_plotting = False
except ImportError as e:
    print("matplotlib is not installed. Plotting is disabled. Matplotlib can be installed via\n pip install matplotlib\n")
    disable_plotting = True

try:
    from tqdm import tqdm
except ImportError as e:
    print("tqdm is not installed. The program will not show progress bars when computing. Tqdm can be installed via\n pip install tqdm\n")

    # work around
    def tqdm(iterator=None, total=None, leave=None):
        if iterator is not None:
            return iterator
        else:
            class Stub:
                def update(self, *args, **kwargs):
                    return

                def __enter__(self):
                    return

                def __exit__(self, exc_type, exc_val, exc_tb):
                    return
            return Stub()



_INFINITY = 10000000000
_BRILL_NOETHER = 1000
_HILBERT_CHOW = 1001
_LI_GIESEKER_UHLENBECK = 1002

MAX_DEPTH_DIVISORIAL = 2  # Maximal depth of search for divisorial
MAX_DEPTH_SPHERICAL_FLOPS = 2  # Maximal depth of search for spherical flops
MAX_DEPTH_SUM_FLOPS = 1  # Maximal depth of search for sum flops


class ModuliSpace:
    """
    Class that implements an instance of a moduli space M(r,D,c)
    for a K3 and all functionality to compute the walls
    """
    def __init__(self, k3: K3Surface, mukai_vec: mVec, printing: bool = True) -> None:
        """
        :param k3: K3 surface
        :param mukai_vec: Mukai vector for the modulispace, i.e. (r,D,c) for M(r,D,c).
        :param printing: Prints steps of computation (True)
        """
        self._k3 = k3
        self._mukai_vec = mukai_vec
        self._printing = printing
        # print(f"DEPTH: {MAX_DEPTH_DIVISORIAL}")


    def get_dim(self) -> int:
        """
        :return: Dimension of the Moduli space, i.e. v^2+2.
        """
        return self.get_mukai() * self.get_mukai() + 2

    def get_mukai(self) -> mVec:
        return self._mukai_vec

    def get_k3(self) -> K3Surface:
        return self._k3

    @staticmethod
    def plot_graph(list_of_num_walls: list) -> None:
        """
        Plot walls in 2D alpha-beta-plane, Walls are the circles or lines with v(vec1) = v(vec2).

        :param list_of_num_walls: List of tuples of vectors or (vec1, vec2, colour, name= vec1)
        :return: None
        """
        if disable_plotting:
            print("Plotting disabled as matplotlib is not installed. Install it with \n pip install matplotlib")
            return

        fig, ax = plt.subplots()
        plt.xlabel("Beta")
        plt.ylabel("Alpha")

        max_size = 0
        min_size = 0

        max_h = 0

        for wall in list_of_num_walls:
            vec1 = wall[0]
            vec2 = wall[1]
            color = 'black'
            if len(wall) >= 3:
                color = wall[2]
            if len(wall) >= 4:
                pass

            (pre_q, pre_l, pre_c) = vec1.get_quadratic_wall_equation(vec2)

            if pre_q == 0:
                if pre_l == 0:
                    # there is no x term in the equation, i.e. only c == 0
                    continue
                x = (pre_c/pre_l).evalf()
                line = plt.axvline(x, color=color)
                print("Vertical: beta = " + str(pre_c/pre_l))
                ax.add_patch(line)

                if x > max_size:
                    max_size = x
                if x < min_size:
                    min_size = x
                continue

            center, rad_squared = ModuliSpace._get_center_radius_squared(pre_q, pre_l, pre_c)
            if center is Rational:
                center = center.evalf()
            if rad_squared is Rational:
                rad_squared = rad_squared.evalf()

            if rad_squared >= 0:
                sqroot = math.sqrt(rad_squared)
                print("Circle: ", center, sqroot)
                circle = plt.Circle((center, 0), sqroot, color=color, fill=False)
                if sqroot > max_h:
                    max_h = sqroot
                if center-sqroot < min_size:
                    min_size = center-sqroot
                if center+sqroot > max_size:
                    max_size = center-sqroot
                ax.add_patch(circle)

        if max_h == 0:
            max_h = 1.
        if min_size is Rational:
            min_size = min_size.evalf()
        if max_size is Rational:
            max_size = max_size.evalf()
        if max_h is Rational:
            print("max_h cannot be Rational.")
            # max_h = max_h.evalf()
        min_size = float(min_size)
        max_size = float(max_size)
        max_h = float(max_h)
        plt.xlim([min_size-0.5, max_size+0.5])
        plt.ylim([0.0, max_h+0.5])

        plt.show()

    @staticmethod
    def _get_center_radius_squared(pre_q: Rational,
                                   pre_l: Rational,
                                   pre_c: Rational) -> tuple[Rational, Rational]:
        center = -pre_l/(2*pre_q)
        rad_sq = pre_l*pre_l/(4*pre_q*pre_q)-pre_c/pre_q
        return center, rad_sq

    def compute_walls(self, ample: mVec, spherical: bool = True,
                      spherical_codims: Optional[list[int]] = None, sum_flops: bool = True) -> list[Wall]:
        """
        Compute the walls. Assumes that boundary of the effective cone is rational.

        :param ample: Specify an ample class for the moduli space.
        :param spherical: Computes spherical walls (True)
        :param spherical_codims: List of codimensions of shperical flops to compute. If None all codimensions.
        :param sum_flops: Compute sum flops (True)
        :return: list of walls.
        """
        # assumes that effective boundary is rational
        walls = self.get_movable_boundary(ample)
        walls = list(walls)

        if walls[0] is None or walls[1] is None:
            print("Did not find divisorial contractions or Jacobians. Either no walls exist or MAX_DEPTH is too small.")
            raise NotImplementedError("Finding flops not implemented yet if boundary divisors are irrational.")

        mov_1 = walls[0].get_movable_divisor()
        mov_2 = walls[1].get_movable_divisor()

        if spherical:
            if spherical_codims is None:
                for i in range(Rational(self.get_dim(), 1) / 2 - 1):
                    walls.extend(self.get_flops_spherical(mov_1, mov_2, i + 2))
            else:
                for i in spherical_codims:
                    walls.extend(self.get_flops_spherical(mov_1, mov_2, i + 2))
        if sum_flops:
            walls.extend(self.get_flops_sum(mov_1, mov_2))
        return walls

    def order_walls_basis(self, walls: list[Wall], basis_1: mVec, basis_2: mVec) -> list[Wall]:
        """
        Order walls corresponding to the basis vectors.

        :return: Sorted list of walls.
        """
        walls.sort(key=lambda x: self._get_wall_thres(x, basis_1, basis_2))
        return walls

    def _get_wall_thres(self, wall: Wall, basis_1: mVec, basis_2: mVec) -> Rational:
        """
        For movable = xB1 + yB2 returns x/y. This is used in order to get an ordering
        of vectors w.r.t. to a chosen basis of the lattice, see order_walls_basis().
        """
        mov = wall.get_movable_divisor()
        return self._get_div_thres(mov, basis_1, basis_2)

    def _get_div_thres(self, div: mVec, basis_1: mVec, basis_2: mVec) -> Rational:
        """
        For div = xB1 + yB2 returns x/y. This is used in order to get an ordering
        of vectors w.r.t. to a chosen basis of the lattice.
        """
        x, y = self._get_factors(div, basis_1, basis_2)
        if y == 0:
            return Rational(_INFINITY, 1)
        return x/y

    def print_movable_divisors(self, walls: list[Wall], basis_1: mVec, basis_2: mVec, factor: float = 1, rational: bool = True,
                               name_basis_1: str = "B_1", name_basis_2: str = "B_2", print_contracted: bool = False) -> None:
        """
        Prints movable divisor of walls in the specified basis.
        Also prints the contracted ones if print_contracted ==True.
        Factor is an artificial multiplicative factor shown in front of the first basis_1.
        """
        walls = self.order_walls_basis(walls, basis_1, basis_2)
        for wall in walls:
            thres = self._get_wall_thres(wall, basis_1, basis_2)
            if not rational:
                thres = thres.evalf()
            basis_str = " "+name_basis_1+"+"+name_basis_2+"\t "
            string = str(thres*factor) + basis_str + wall.get_type_str()
            if thres == _INFINITY:
                basis_str = " "+name_basis_1 + "\t\t "
                string = "1" + basis_str + wall.get_type_str()
                if not rational:
                    string = "1.0000000" + basis_str + wall.get_type_str()
            if wall.get_type() == Wall.FLOP_SPHERE or wall.get_type() == Wall.FL_SUM:
                string += " of codimension "+str(wall.get_codim())
            elif wall.get_type() in [Wall.BN, Wall.LGU, Wall.HC] and print_contracted:
                div = wall.get_contracted_divisor()
                thres_div = self._get_div_thres(div, basis_1, basis_2)
                string_div = ""
                if not rational:
                    thres_div = thres.evalf()
                if thres_div != _INFINITY:
                    string_div = f" contracting {thres_div} {name_basis_1} + {name_basis_2}\t"
                else:
                    thres_div = "1"
                    if not rational:
                        thres_div = "1.0000000"
                    string_div = f" contracting {thres_div} {name_basis_1} + 0 {name_basis_2}"
                string += string_div

            print(string)
        return

    def compute_walls_printing(self, ample: mVec, basis_1: mVec, basis_2: mVec,
                               rational: bool = True, spherical: bool = True,
                               spherical_codims: Optional[list[int]] = None,
                               sum_flops: bool = True, factor: float = 1,
                               name_1: str = "B_1", name_2: str = "B_2") -> list[Wall]:
        """
        Computes walls and prints walls in the specified basis.
        Factor is an artificial multiplicative factor shown in front of the first basis_1.

        :param ample: Ample divisor for moduli space
        :param basis_1:
        :param basis_2:
        :param rational: Display numbers as fractions
        :param spherical: Compute spherical flops
        :param spherical_codims: If None compute all sperical flops, otherwise List of codimensions to compute.
        :param sum_flops: Compute sum flops (True)
        :param factor: Mupltiply the first factor in the basis representation.
        :param name_1: Names of the basis_1
        :param name_2: Names of the basis_2
        :return: List of walls.
        """
        walls = self.compute_walls(ample, spherical=spherical, spherical_codims=spherical_codims, sum_flops=sum_flops)
        self.print_movable_divisors(walls, basis_1, basis_2, factor=factor, rational=rational,
                                    name_basis_1=name_1, name_basis_2=name_2)
        return walls

    def _get_equation_vector_given_self_and_mukai_intersection(self, self_intersection: int,
                                                               intersection: int) -> QEquation:
        """
        Returns sympy equation in coordinates x,y for w^2 = self_intersection, w*v = intersection.
        Solutions are for a vector w = (y,x, ...).
        """
        muk = self.get_mukai()
        k3 = self.get_k3()
        a = muk.get_rank() * k3.get_deg()
        b = 2 * muk.get_ch2()
        c = -2 * muk.get_divisor() * k3.get_deg()
        d = 2 * intersection
        e = -self_intersection * muk.get_rank()
        return QEquation(a, b, c, d, e)

    def _get_isotropic_with_intersection(self, intersection: int) -> QEquation:
        """
        Returns equation for w^2 = 0, w*v = intersection
        Solutions are for a vector w = (y,x, ...).
        """
        return self._get_equation_vector_given_self_and_mukai_intersection(0, intersection)

    def _get_spherical_with_given_intersection(self, intersection: int) -> QEquation:
        """
        Returns equation for w^2 = -2, w*v = intersection
        Solutions are for a vector w = (y,x, ...).
        """
        return self._get_equation_vector_given_self_and_mukai_intersection(-2, intersection)

    def _get_spherical_flop_equation(self, codim: int) -> QEquation:
        """
        Returns equation for w^2 = -2, w*v = codim-1
        Solutions are for a vector w = (y,x, ...).
        """
        return self._get_spherical_with_given_intersection(codim - 1)

    def _check_bounds(self, vec: mVec, ample: mVec, ample_orthogonal: mVec, pos_thres: float, neg_thres: float) -> bool:
        """Check whether vector has thres against ample in between pos_thres and neg_thres."""
        thres = vec.get_thres(ample, ample_orthogonal)
        if pos_thres > thres > neg_thres:
            return True
        else:
            return False

    def _get_factors(self, vec: mVec, basis_1: mVec, basis_2: mVec) -> tuple[Rational, Rational]:
        """
        Computes alpha beta in vec = alpha*basis_1+beta*basis_2.

        :return: alpha, beta
        """
        matrix = Matrix(
            [[basis_1.get_rank(), basis_2.get_rank()], [basis_1.get_divisor(), basis_2.get_divisor()],
             [basis_1.get_ch2(), basis_2.get_ch2()]])
        b = Matrix([[vec.get_rank()], [vec.get_divisor()], [vec.get_ch2()]])

        sol, params = matrix.gauss_jordan_solve(b)
        if len(params) > 0:
            print("Boundaries were not independent.")
            raise ValueError
        return sol[0], sol[1]

    def _check_in_cone(self, vec: mVec, bound_1: mVec, bound_2: mVec) -> bool:
        """Checks whether vec is contained in the cone spanned by bound_1, bound_2."""
        x, y = self._get_factors(vec, bound_1, bound_2)
        if x > 0 and y > 0:
            return True
        return False

    def vec_induces_divisorial(self, vec: mVec) -> bool:
        """
        Checks whether lattice generated by vec and mukai vector of moduli space
        induces a divisorial contraction.
        """
        muk = self.get_mukai()
        vec_2 = vec*vec
        muk_2 = muk*muk
        vec_m = vec*muk
        det = vec_2*muk_2 - vec_m*vec_m
        if not det < 0:  # Check if Lattice is hyperbolic
            return False

        if is_square(-det):
            # check if lattice induces an LGU or HC contraction
            # How: Computes isotropic final_vec which has final_vec.muk = 2
            # If -det is not a square the vector can not possibly be rational
            x = Rational(2, 1)/integer_nthroot(-det, 2)[0]
            y = (2 - x * vec_m) / muk_2
            final_vec = x*vec + y*muk
            if final_vec.is_integral():
                return True
            y = (2 + x * vec_m) / muk_2
            final_vec = -x * vec + y * muk
            if final_vec.is_integral():
                return True

        # Check if BN:
        # How: vector must be orthogonal to m (Mukai).
        # ->v' =  v - vm/m^2 m.
        # (beta v')^2 = -2
        # Check if beta is square and then comput if beta v' is integral.
        sq = -2*muk_2/det
        sq = simplify(sq)
        n, d = sq.as_numer_denom()

        if is_square(n) and is_square(d):
            x = Rational(integer_nthroot(n, 2)[0], integer_nthroot(d, 2)[0])
            y = Rational(-x * vec_m, 1) / muk_2
            final_vec = x * vec + y * muk
            if final_vec.is_integral():
                return True
            y = x * vec_m / muk_2
            final_vec = -x * vec + y * muk
            if final_vec.is_integral():
                return True

        return False

    def get_flops_sum(self, mov_1: mVec, mov_2: mVec) -> list[Wall]:
        """Computes the flopping contractions that arise as v = v_1+v_2,
        where v_1, v_2 are positive in the sense of Bayer--Macrí,
        i.e. that v^2 >0 and in the component of v.mukai > 0,
        such that furthermore the corresponding movable divisor is inside the cone spanned by mov_1, mov_2.

        :return: List of walls.
        """
        if self._printing:
            print("Compute Flops v = a+b.")
        flops = []
        muk = self.get_mukai()
        ample = mov_1+mov_2  # not necessarily ample, but is in the inside of the movable cone

        muk_2 = Rational(muk*muk, 1)

        # with tqdm(total=int((muk_2/2) * (muk_2/2-1))) as pbar:
            # Go through all possibilities of v = vec+vec2 such that vec^2 = 2k, vec*muk-1 = m.
        for m in tqdm(range(muk_2/2)):    # v*mukai-1
            for k in tqdm(range(muk_2/2-1), leave=False):  # n = v^2 up to muk^2-4 (as otherwise v^2 > muk^2).
                n = 2*k  # as n = v^2 is even

                eq = self._get_equation_vector_given_self_and_mukai_intersection(n, m + 1)
                solutions = eq.get_solutions_in_range(MAX_DEPTH_SUM_FLOPS)
                for a_final, b_final in solutions:
                    vec = self._get_integral_vec_with_given_intersection(m + 1, a_final, b_final)
                    if vec is None:
                        continue
                    vec_2 = muk-vec
                    if vec_2*vec_2 < 0:
                        # check whether both classes are positive in the sense of Bayer--Macrí.
                        # vec*vec > 0 and vec*muk > 0 is ensured by construction.
                        continue
                    if self.vec_induces_divisorial(vec):
                        continue

                    codim = vec*vec_2
                    if not codim.is_integer:
                        raise Exception(f"Something went wrong: codim {codim} should be an Integer")
                    codim = int(codim)

                    vec_final = vec.get_orthogonal(muk, ample)
                    if self._check_in_cone(vec_final, mov_1, mov_2):
                        wall = Wall(vec, Wall.FL_SUM, vec_final, codim=codim)
                        wall.make_primitive()
                        if wall not in flops:
                            flops.append(wall)
                # pbar.update(1)
        return flops

    def get_flops_spherical(self, mov_1: mVec, mov_2: mVec, codim: int) -> list[Wall]:
        """
        returns spherical flops that are in the chamber between mov_1, mov_2
        I.e. finds classes such s that s^2 = -2 and s.v = codim-1.
        Then checks whether orthogonal of s is in between the two movable divisors.
        """
        if self._printing:
            print("Compute Flops v^2 = -2 with codim " + str(codim))

        flops = []
        flop_walls = []
        muk = self.get_mukai()
        ample = mov_1+mov_2  # not necessarily ample but movable, which suffices

        # Spherical flops
        # s^2 = -2 , s.v = codim -1
        eq = self._get_spherical_with_given_intersection(codim - 1).get_solutions_in_range(MAX_DEPTH_SPHERICAL_FLOPS, tqdm_=tqdm)
        if len(eq) == 0:  # no solution
            return flops

        for x_final, y_final in tqdm(eq):
            vec = self._get_integral_vec_with_given_intersection(codim - 1, x_final, y_final)  # spherical class
            if vec is None:
                continue

            lattice = Lattice(self.get_mukai(), vec)
            if lattice.is_hyperbolic() is False or self.vec_induces_divisorial(vec):
                # hyperbolic should follow from intersection numbers
                continue

            vec_final = vec.get_orthogonal(muk, ample).get_primitive()  # line bundle

            if vec_final in flops:
                continue

            if self._check_in_cone(vec_final, mov_1, mov_2):
                flops.append(vec_final)
                wall = Wall(vec, Wall.FLOP_SPHERE, vec_final, codim=codim)
                flop_walls.append(wall)
            else:
                continue
        return flop_walls

    def get_movable_boundary(self, ample: mVec) -> tuple[Optional[Wall], Optional[Wall]]:
        """Returns Walls that are the boundary of the movable cone if they exist.
        A Lagrangian fibration is also regarded as wall, contrary to Bayer--Macrí.

        Otherwise, they are irrational.

        That is, it computes Brill-Noether, Hilbert-Chow and Li-Gieseker-Uhlenbeck contractions.
        That is, it finds classes s with
            - s^2 = -2 and s.v = 0 (BN)
            - s^2 = 0 and s.v = 1 (HC)
            - s^2 = 0 and s.v = 2 (LGU)

        Hilbert-Chow and Li-Gieseker-Uhlenbeck are computed in one go, by checking if s is divisible by 2.

        The returned contracted divisors are a multiple of the actual one which is indivisible.
        """
        if self._printing:
            print("Compute movable cone.")

        boundary = []
        boundary_walls = []

        eq = {}
        orthogonal = ample.get_orthogonal(self.get_mukai())

        intersections = {_BRILL_NOETHER: 0, _HILBERT_CHOW: 1, _LI_GIESEKER_UHLENBECK: 2}
        wall_types = {_BRILL_NOETHER: Wall.BN, _HILBERT_CHOW: Wall.HC, _LI_GIESEKER_UHLENBECK: Wall.LGU}

        for flopping_type in [_BRILL_NOETHER, _HILBERT_CHOW, _LI_GIESEKER_UHLENBECK]:
            intersection = intersections[flopping_type]

            if flopping_type == _BRILL_NOETHER:
                if self._printing:
                    print("Compute Brill-Noether contractions")
                eq = self._get_spherical_with_given_intersection(0)

            elif flopping_type == _HILBERT_CHOW:
                # This case is obsolete now, as HC can be computed in LGU case by %2 = 0.
                continue

            elif flopping_type == _LI_GIESEKER_UHLENBECK:
                if self._printing:
                    print("Compute Hilbert-Chow and Li-Gieseker-Uhlenbeck contractions")
                eq = self._get_isotropic_with_intersection(2)

            list_solutions = eq.get_solutions_in_range(MAX_DEPTH_DIVISORIAL, tqdm_=tqdm)

            for (a_final, b_final) in tqdm(list_solutions):
                vec = self._get_integral_vec_with_given_intersection(intersection, a_final, b_final)

                if vec is None:
                    continue

                if flopping_type == _BRILL_NOETHER:
                    lattice = Lattice(self.get_mukai(), vec)
                    if lattice.contains_isotropic_element():
                        # Then it is a HC or LGU
                        continue

                vec_final = vec
                wall_type = wall_types[flopping_type]  # only changed when LGU is actually HC

                if flopping_type == _LI_GIESEKER_UHLENBECK:
                    # check that it is really an LGU and not HC
                    if vec.get_rank() % 2 == 0 and vec.get_divisor() % 2 == 0 and vec.get_ch2() % 2 == 0:
                        # if HC then take primitive of vec_final as this will have vec_final.muk = 1
                        vec = vec.get_primitive()
                        wall_type = Wall.HC  # i.e. it is HC

                # if not LGU or HC done
                if flopping_type in [_LI_GIESEKER_UHLENBECK, _HILBERT_CHOW]:
                    vec_final = (self.get_mukai() * self.get_mukai()) * vec + (-(self.get_mukai() * vec)) * self.get_mukai()

                if vec_final*ample <= 0:
                    vec_final = -vec_final
                vec_final = vec_final.get_primitive()  # the contracted divisor

                if vec_final in boundary:
                    # Check if there is already a wall that trumps the current wall:
                    # i.e. HC > LGU > BN
                    # and only one of these is permitted
                    walls_with_given_vec = [w.get_type() for w in boundary_walls if w.get_contracted_divisor() == vec_final]
                    if Wall.HC in walls_with_given_vec:
                        continue
                    if wall_type == Wall.BN:
                        continue
                    if wall_type in walls_with_given_vec:
                        continue

                    boundary_walls = [w for w in boundary_walls if w.get_contracted_divisor() != vec_final]
                    boundary = [b for b in boundary if b != vec_final]

                movable_divisor = vec_final.get_orthogonal(self.get_mukai(), ample)
                wall = Wall(vec, wall_type, movable_divisor, contracted_divisor=vec_final)
                boundary_walls.append(wall)
                boundary.append(vec_final)


        # get the one with the smallest distance

        pos_final = None
        neg_final = None
        thres_pos = None
        thres_neg = None
        for wall in boundary_walls:
            vec = wall.get_movable_divisor()

            thres = vec.get_thres(ample, orthogonal)

            if thres > 0:
                if pos_final is None:
                    pos_final = wall
                    thres_pos = thres
                elif thres < thres_pos:   # get smallest thres
                    pos_final = wall
                    thres_pos = thres
            elif thres < 0:
                if neg_final is None:
                    neg_final = wall
                    thres_neg = thres
                elif thres > thres_neg:
                    neg_final = wall
                    thres_neg = thres

        # If one of the walls is not found yet, look for Lagrangians, i.e. divisors with D^2 = 0:
        if neg_final is None or pos_final is None:
            jacobians = []
            eq = self._get_isotropic_with_intersection(0)
            for t in eq.get_equation():
                a = t[0]
                b = t[1]
                variables = a.free_symbols
                if len(variables) == 0:
                    variables = b.free_symbols
                if len(variables) > 0:
                    var = variables.pop()
                    a = simplify(a.subs(var, 1))
                    b = simplify(b.subs(var, 1))
                if (a, b) == (0, 0):
                    continue
                vec_final = self._get_vec_with_given_intersection(0, a, b)
                # as isotropic, integrality is not needed

                if vec_final*ample <= 0:
                    vec_final = -vec_final
                if vec_final not in jacobians:
                    jacobians.append(vec_final)
            for vec in jacobians:
                thres = vec.get_thres(ample, orthogonal)
                wall = Wall(vec.get_primitive(), Wall.LAGRANGIAN, vec.get_primitive())
                if thres < 0 and neg_final is None:
                    neg_final = wall
                if thres > 0 and pos_final is None:
                    pos_final = wall

        if pos_final is not None:
            pos_final.make_primitive()
        if neg_final is not None:
            neg_final.make_primitive()

        return pos_final, neg_final

    def _get_integral_vec_with_given_intersection(self, intersection: int,
                                                  x: Union[Rational, int],
                                                  y: Union[Rational,  int]) -> Optional[mVec]:
        """
        Returns vector (y,x, ...) with given Intersection against Mukai vector
        Returns None, if vector (y,x, ....) is not integral
        """

        c = -intersection + x * self.get_k3().get_deg() * self.get_mukai().get_divisor() - y * self.get_mukai().get_ch2()
        if self.get_mukai().get_rank() == 0:
            if self.get_mukai().get_divisor() * x * self.get_k3().get_deg() - y * self.get_mukai().get_ch2() != intersection:
                return None
            else:
                return mVec(self.get_k3(), y, x, Rational(0, 1))
        if c % self.get_mukai().get_rank() != 0:
            return None
        return mVec(self.get_k3(), y, x, c / self.get_mukai().get_rank())

    def _get_vec_with_given_intersection(self, intersection: int,
                                         x: Union[Rational, int],
                                         y: Union[Rational, int]) -> mVec:
        """Returns vector (y,x, ...) with given Intersection against Mukai vector"""
        return mVec(self.get_k3(), y, x, (
                -intersection + x * self.get_k3().get_deg() * self.get_mukai().get_divisor() - y * self.get_mukai().get_ch2()) / self.get_mukai().get_rank())

    def _get_orthogonal_to_mukai(self, vec: mVec) -> mVec:
        """Returns vector which is orthogonal to both vec and mukai vector of Moduli space"""
        return self.get_mukai().get_orthogonal(vec)

    def __str__(self) -> str:
        muk = self.get_mukai()
        string = "M" + str(muk)
        if muk.get_rank() == 1 and muk.get_divisor() == 0:
            string = "Hilb^[" + str(-muk.get_ch2() + 1) + "]"
        return string


def write_walls_file(file: TextIO, modulispace: ModuliSpace,
                     list_walls: list[Wall], ample: Optional[mVec] = None) -> TextIO:
    """
    Writes list_walls into the file object.

    :param file: Opened file object
    :param modulispace: The modulispace the walls correspond to
    :param list_walls: The walls
    :param ample: The ample divisor for the walls that was specified during computation
    :return: the file object
    """
    file.write(str(modulispace.get_k3().get_deg()) + "\n")
    file.write(str(modulispace.get_mukai()) + "\n")
    if ample:
        file.write("#"+str(ample))
        file.write("\n")
    for wall in list_walls:
        wall.make_primitive()
        file.write(str(wall.get_vec()))
        file.write(";")
        file.write(str(wall.get_type()))
        file.write(";")
        file.write(str(wall.get_movable_divisor()))
        file.write(";")
        file.write(str(wall.get_contracted_divisor()))
        file.write(";")
        file.write(str(wall.get_codim()))
        file.write("\n")
    return file


def _get_vec_from_string(k3: K3Surface, string: str) -> mVec:
    """
    Returns the mVec for a String of type '(r,k,c)'

    :param k3: K3 surface
    :param string: String of the form '(r,k,c)'
    :return: mVec(k3, r, k, c)
    """
    vec_str = string[1:-1]
    nums = vec_str.split(",")
    return mVec(k3, Rational(int(nums[0]), 1), Rational(int(nums[1]), 1), Rational(int(nums[2]), 1))


def read_walls_file(file: TextIO, return_ample: bool = False):
    """
    Reads walls in the opened file object.
    Does not check for correctness of file.

    :param file: file object opened for reading
    :param return_ample:
    :return: modulispace, list of walls, ample (if return_ample==True)
    """
    list_walls = []
    line1 = file.readline()
    line1 = line1.replace("\n", "")
    k3 = K3Surface(int(line1))
    line2 = file.readline()
    line2 = line2.replace("\n", "")
    muk = _get_vec_from_string(k3, line2)
    moduli = ModuliSpace(k3, muk)
    ample = None

    for line in file:
        line = line.replace("\n", "")
        if line[0] == '#':
            line = line.replace("#", "")
            ample = _get_vec_from_string(k3, line)
            continue
        vecs = line.split(";")
        special_vec = _get_vec_from_string(k3, vecs[0])
        wall_type = int(vecs[1])
        divisor_vec = _get_vec_from_string(k3, vecs[2])
        if vecs[3] != "None":
            contracted_vec = _get_vec_from_string(k3, vecs[3])
        else:
            contracted_vec = None
        if vecs[4] != "None":
            codim = int(vecs[4])
        else:
            codim = None
        wall = Wall(special_vec, wall_type, divisor_vec, contracted_vec, codim)
        list_walls.append(wall)

    if return_ample:
        if ample is not None:
            return moduli, list_walls, ample
        else:
            raise Exception(f"File '{file.name}' does not contain an ample vector!")
    return moduli, list_walls
