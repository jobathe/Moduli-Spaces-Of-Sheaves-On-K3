"""
This file comprises the functionality of a Mukai vector with its intersection pairing.

Mukai = (a,k,b) in Z + Pic X + Z.
(a,k,b).(c,k',d) = kk'-ad-bc

It also presents the functionality of:
    - Computing the Z function for Bridgeland stability conditions in the alpha-beta plane.
    - Computing the quadratic equation for a wall, i.e. the equation coming from v(self) = v(vec)
      in the alpha-beta-plane.
    - Computing orthogonals.
    - Computes whether the vector is effective in the sense of Bayer--Macrí.
    - Computation of "threshold" which induces an ordering of vectors w.r.t. some fixed basis in a lattice.
"""
from __future__ import annotations

from sympy.matrices import Matrix
from sympy.core.numbers import Rational
from sympy.core.numbers import igcd
from sympy.core.numbers import ilcm

from modulispace.k3surface import K3Surface

from typing import Union, Optional


class mVec:
    def __init__(self, k3: K3Surface, rk: Union[int, Rational],
                 divisor: Union[int, Rational], ch2: Union[int, Rational]):
        self.k3 = k3

        if not isinstance(rk, Rational):
            rk = Rational(rk, 1)
        if not isinstance(divisor, Rational):
            divisor = Rational(divisor, 1)
        if not isinstance(ch2, Rational):
            ch2 = Rational(ch2, 1)
        self.rk = rk
        self.div = divisor
        self.ch2 = ch2

    def get_k3(self) -> K3Surface:
        return self.k3

    def get_rank(self) -> Union[int, Rational]:
        return self.rk

    def get_divisor(self) -> Union[int, Rational]:
        return self.div

    def get_ch2(self) -> Union[int, Rational]:
        return self.ch2

    def is_effective(self, mukai: mVec, alpha: Rational, beta: Rational) -> bool:
        """
        is this vector effective w.r.t. to mukai vector, i.e. Re(Z(self)/Z(mukai)) > 0
        """
        r_self, i_self = self.get_z(alpha, beta)
        r_muk, i_muk = mukai.get_z(alpha, beta)

        const = r_self*r_muk+i_self*i_muk

        if const > 0:
            return True
        return False

    def get_z(self, alpha: Rational, beta: Rational) -> tuple[Rational, Rational]:
        """
        Computes the complex Z function from Bayer--Macrí for alpha-beta in the alpha-beta-plane.
        :param alpha:
        :param beta:
        :return: r,i (where Re(Z) = r, Im(Z)=i)
        """
        ch0 = self.get_rank()
        ch1 = self.get_divisor()
        ch2 = self.get_ch2() - ch0
        deg = self.get_k3().get_deg()

        r = -ch2 + ch1*beta*deg + (-beta*beta*ch0 + alpha*alpha*ch0)/2
        i = alpha*ch1*deg - alpha*beta*deg

        return r, i

    def get_quadratic_wall_equation(self, vec: mVec) -> tuple[Rational, Rational, Rational]:
        """
        Returns factors for a_2x^2+a_1x+a_0 = 0, which is the equation of v(self) = v(vec)
        :param vec:
        :return: a_2, a_1, a_0
        """
        ch0s = self.get_rank()
        ch1s = self.get_divisor()
        ch2s = self.get_ch2() - ch0s

        ch0v = vec.get_rank()
        ch1v = vec.get_divisor()
        ch2v = vec.get_ch2() - ch0v

        deg = self.get_k3().get_deg()

        pre_square = -deg*(ch0s*ch1v-ch1s*ch0v)/2
        pre_lin = ch2v*ch0s-ch0v*ch2s
        pre_const = ch1v*ch2s-ch1s*ch2v

        return pre_square, pre_lin, pre_const

    def vec_prod(self, other: mVec) -> Rational:
        return self.get_divisor() * other.get_divisor() + self.get_ch2() * other.get_ch2() + other.get_rank() * self.get_rank()

    def project_to_vec(self, other: mVec) -> Rational:
        return (self*other)*other

    def get_thres(self, ample: mVec, orthogonal: mVec) -> Optional[Rational]:
        """
        computes the threshold against ample and orthogonal,
        i.e. self = x ample + y orthogonal
        If ample, orthogonal are colinear returns None

        :returns y/x
        """

        matrix = Matrix([[ample.get_rank(), orthogonal.get_rank()],
                         [ample.get_divisor(), orthogonal.get_divisor()],
                         [ample.get_ch2(), orthogonal.get_ch2()]])
        b = Matrix([[self.get_rank()], [self.get_divisor()], [self.get_ch2()]])

        sol, params = matrix.gauss_jordan_solve(b)

        if len(params) > 1:
            return None

        if len(params) == 0:
            x = sol[0]
            y = sol[1]
        else:
            var = params[0]
            specific_sol = sol.replace(var, 1)
            x = specific_sol[0]
            y = specific_sol[1]

        return y/x

    def get_vec_orthogonal(self, mukai: mVec) -> Optional[mVec]:
        """
        get vector orthogonal to mukai in intersection pairing and orthogonal to self in euclidean
        None, if mukai and self are colinear.
        """

        d = self.get_k3().get_deg()

        matrix = Matrix(
            [[self.get_rank(), self.get_divisor(), self.get_ch2()],
             [-mukai.get_ch2(), d * mukai.get_divisor(), -mukai.get_rank()]])
        b = Matrix([[0], [0]])

        sol, params = matrix.gauss_jordan_solve(b)

        if len(params) > 1:
            return None

        var = params[0]

        specific_sol = sol.replace(var, 1)

        return mVec(self.get_k3(), specific_sol[0], specific_sol[1], specific_sol[2])

    def get_primitive(self) -> mVec:
        """Returns primitive form of vector, i.e. the smallest integral vector that is a multiple of self"""
        r = self.get_rank()
        d = self.get_divisor()
        c = self.get_ch2()

        r_n, r_d = r.as_numer_denom()
        d_n, d_d = d.as_numer_denom()
        c_n, c_d = c.as_numer_denom()

        lcm = ilcm(r_d, d_d, c_d)
        r = lcm*r
        d = lcm*d
        c = lcm*c
        gcd = igcd(r, d, c)
        return mVec(self.get_k3(), r / gcd, d / gcd, c / gcd)

    def __mul__(self, vec2: mVec) -> Rational:
        return self.__rmul__(vec2)

    def __rmul__(self, d: Union[Rational, int, mVec]) -> Union[Rational, mVec]:
        """
        Returns vector multiplication if d is a vector, and scalar multiplication otherwise.
        """
        if isinstance(d, mVec):
            return self.k3.get_deg() * self.get_divisor() * d.get_divisor() - self.get_rank() * d.get_ch2() - d.get_rank() * self.get_ch2()
        else:
            return mVec(self.get_k3(), self.get_rank() * d,
                        self.get_divisor() * d, self.get_ch2() * d)

    def __add__(self, vec2: mVec) -> mVec:
        return mVec(self.k3, self.get_rank() + vec2.get_rank(),
                    self.get_divisor() + vec2.get_divisor(),
                    self.get_ch2() + vec2.get_ch2())

    def __sub__(self, vec2: mVec) -> mVec:
        return mVec(self.k3, self.get_rank() - vec2.get_rank(),
                    self.get_divisor() - vec2.get_divisor(),
                    self.get_ch2() - vec2.get_ch2())

    def __str__(self) -> str:
        return str((self.get_rank(), self.get_divisor(), self.get_ch2()))

    def __neg__(self) -> mVec:
        return mVec(self.get_k3(), -self.get_rank(), -self.get_divisor(), -self.get_ch2())

    def __eq__(self, other: mVec) -> bool:
        return self.get_rank() == other.get_rank() and self.get_divisor() == other.get_divisor() and self.get_ch2() == other.get_ch2()

    def __repr__(self) -> str:
        return self.__str__()

    def is_integral(self) -> bool:
        return self.get_rank().is_integer and self.get_divisor().is_integer and self.get_ch2().is_integer

    def get_orthogonal(self, vec2: mVec, ample: Optional[mVec] = None):
        """
        Returns orthogonal to self and vec2. If ample is not none, the result satisfies ample*result > 0.
        """
        d = self.get_k3().get_deg()

        matrix = Matrix(
            [[-self.get_ch2(), d * self.get_divisor(), -self.get_rank()],
             [-vec2.get_ch2(), d * vec2.get_divisor(), -vec2.get_rank()]])
        b = Matrix([[0], [0]])

        sol, params = matrix.gauss_jordan_solve(b)

        if len(params) > 1:
            return None

        var = params[0]

        specific_sol = sol.replace(var, 1)

        vec_final = mVec(self.get_k3(), specific_sol[0], specific_sol[1], specific_sol[2])
        if ample is not None:
            if vec_final*ample < 0:
                vec_final = -vec_final

        return vec_final
