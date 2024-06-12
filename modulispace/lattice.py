from sympy import symbols
from sympy.solvers.diophantine import diophantine
from modulispace.mukaivector import mVec


class Lattice:
    """
    Represents a lattice spanned by two vectors.
    """
    def __init__(self, vec1: mVec, vec2: mVec) -> None:
        self.vec1 = vec1
        self.vec2 = vec2

    def get_vec1(self) -> mVec:
        return self.vec1

    def get_vec2(self) -> mVec:
        return self.vec2

    def get_rank(self) -> int:
        """
        Returns the rank of the lattice spanned by the two vectors.
        """
        r1 = self.get_vec1().get_rank()
        d1 = self.get_vec1().get_divisor()
        c1 = self.get_vec1().get_ch2()

        r2 = self.get_vec2().get_rank()
        d2 = self.get_vec2().get_divisor()
        c2 = self.get_vec2().get_ch2()

        det1 = r1*d2-r2*d1
        det2 = r1*c2-r2*c1
        det3 = d1*c2-d2*c1

        if det1 == 0 and det2 == 0 and det3 == 0:
            return 1
        return 2

    def is_negative_definite(self) -> bool:
        """
        :return: True if the lattice is negative definite, else false
        """
        vec1 = self.get_vec1()
        vec2 = self.get_vec2()
        v1s = vec1*vec1
        if v1s >= 0:
            return False
        det = v1s*(vec2*vec2) - (vec1*vec2)*(vec1*vec2)
        if det <= 0:
            return False
        return True

    def is_hyperbolic(self) -> bool:
        """
        :return: True if the lattice L is hyperbolic, i.e. if det L < 0.
        """
        vec1 = self.get_vec1()
        vec2 = self.get_vec2()
        v1s = vec1 * vec1
        det = v1s * (vec2 * vec2) - (vec1 * vec2) * (vec1 * vec2)
        if det >= 0 or self.get_rank() < 2:
            return False
        return True

    def contains_isotropic_element(self) -> bool:
        """
        :return: True if the lattice contains an element with v^2 = 0.
        """
        v1s = self.get_vec1() * self.get_vec1()
        v2s = self.get_vec2() * self.get_vec2()
        v1v2 = self.get_vec1() * self.get_vec2()
        x, y = symbols("x, y", integer=True)
        eq = diophantine(x**2*v1s + 2*v1v2*x*y + y**2*v2s)
        if len(eq) > 1:
            return True
        # eq has at least trivial solution, hence only one solution is now allowed

        result_tuple = eq.pop()
        if len(result_tuple) <= 1:
            # I.e. factors of x and/or y are zero
            return True
        a = result_tuple[0]
        b = result_tuple[1]
        if a == 0 and b == 0:
            return False
        return True
