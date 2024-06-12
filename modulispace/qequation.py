from sympy import symbols
from sympy.solvers.diophantine import diophantine
from sympy import simplify
from sympy.core.numbers import Rational
from typing import Union



class QEquation:
    """solves aX^2+bY^2+cYX+dY+e = 0"""
    def __init__(self, a: int, b: int, c: int, d: int, e: int):
        x, y = symbols("x,y", integer=True)
        self.t = symbols("t", integer=True)
        self.eq = diophantine(a * x ** 2 + b * y ** 2 + c * y * x + d * y + e, param=self.t)

    def get_equation(self) -> diophantine:
        return self.eq

    def get_variable(self) -> symbols:
        return self.t

    def get_solutions_in_range(self, range_: Union[iter, int], tqdm_=lambda x: x):
        if isinstance(range_, int):
            range_ = range(-range_ + 1, range_)
        list_ = []
        for a, b in tqdm_(self.eq):
            variables = a.free_symbols
            if len(variables) != 0:
                t_a = variables.pop()
                for i in range_:
                    a_final = simplify(a.subs(t_a, i))
                    b_final = simplify(b.subs(t_a, i))
                    list_.append((a_final, b_final))
            else:
                a_final = simplify(a)
                b_final = simplify(b)
                list_.append((a_final, b_final))
        return list_

    def get_next_solution(self, d: int = 0) -> tuple[Rational, Rational]:
        """Returns solution f(d) of equation"""
        a, b = self.get_equation().pop()

        return simplify(a.subs(self.get_variable(), d)), simplify(b.subs({self.get_variable(): d}))
