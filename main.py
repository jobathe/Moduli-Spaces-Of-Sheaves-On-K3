import os
import modulispace
from modulispace import K3Surface
from modulispace import mVec
from modulispace import ModuliSpace
from modulispace import Wall
from modulispace import read_walls_file, write_walls_file
from sympy.core.numbers import Rational

from typing import Optional


STORAGE_PATH_DIR = os.path.join(os.curdir, "Files")

HELP_STRING_INTERACTIVE = """Hilbert scheme programm. Author: Jonas Baltes. 
Computes and shows walls of a moduli space on a K3 surface with Picard rank 1.

Type one of 
- show|compute|showcompute|plot|rational|showwalls|help|quit
or their first letters 
- s|c|sc|p|r|w|h|q.

Explanation:
- Show: Shows already computed movable divisors and walls (if showwalls is not toggled) of specified degrees and length.
        To show several degrees or length, type the numbers with a space in between or with "-".
        I.e. 2 4 6 8 and 2 - 8 have the same meaning.

- Compute: Computes walls and shows movable divisors of specified degrees and length. Caution: Overwrites already computed files.

- ShowCompute: Shows and (if not already computed) computes the movable cone and walls.

- Plot: Plots walls in the alpha-beta-plane with possible additional user-specified walls

- Rational: Switch between fractions and decimal numbers

- Showwalls: Toggles whether walls shall be shown or just the movable divisors

- Help: Shows this help

- Quit: Quits the programm

"""



def exists_hilb(deg: int, n: int) -> bool:
    """
    Checks if already a file in ./Files/  exists such that it represents the Hilbert scheme of
    length n on a K3 surface of degree deg.

    :param deg: Degree of K3 surface
    :param n: Length of Hilbert scheme Hilb^n
    :return: True if file exists
    """
    name = os.path.join(STORAGE_PATH_DIR, f"Degree{deg}Hilb{n}.walls")
    return os.path.exists(name)


def safe_hilb(hilb: ModuliSpace, walls_list: list[Wall], force_overwrite: bool = False):
    """
    Safes the walls computed for the Hilbert scheme. Returns False if file already exists.

    :param hilb: Hilbert scheme
    :param walls_list: List of walls
    """
    deg = int(hilb.get_k3().get_deg())
    n = int(-hilb.get_mukai().get_ch2() + 1)
    name = os.path.join(STORAGE_PATH_DIR, f"Degree{deg}Hilb{n}.walls")
    print(name)
    if os.path.exists(name) and not force_overwrite:
        return False
    f = open(name, "w")
    write_walls_file(f, hilb, walls_list)
    f.close()
    return True


def compute_hilb(d: int, k: int, spherical: bool = True,
                 sum_flop: bool = True) -> tuple[ModuliSpace, list[Wall]]:
    """
    Computes the wall for the Hilbert scheme of length k for a K3 surface of degree d.
    Does not compute spherical flops if spherical==False
    Does not compute sum flops if sum_flop==False

    :param d: Degree of K3
    :param k: Length of Hilbert scheme
    :param spherical: True if spherical flops are computed
    :param sum_flop: True if sum flops are computed
    :return: (Hilbert scheme, List of walls)
    """

    n = k-1
    k3 = K3Surface(d)
    muk = mVec(k3, Rational(1, 1), Rational(0, 1), Rational(-n, 1))
    amp = mVec(k3, 1, -100, n)
    Hilb = ModuliSpace(k3, muk, printing=True)

    print(str(Hilb) + " of H^2 = " + str(k3.get_deg()))
    basis_1 = mVec(k3, Rational(-1, 1), Rational(0, 1), Rational(-n, 1))
    basis_2 = mVec(k3, Rational(0, 1), Rational(-1, 1), Rational(0, 1))

    walls_list = Hilb.compute_walls_printing(amp, basis_1, basis_2,
                                             name_1="B", name_2="H",
                                             rational=True, sum_flops=sum_flop,
                                             spherical=spherical)

    return Hilb, walls_list


def show_plot(deg: int, n: int, list_of_addtional: Optional[list[tuple[mVec, mVec, str]]] = None) -> None:
    """
    Shows the plot of the alpha-beta-plane for walls for the Hilbert scheme of length n
    of a K3 of degree d.

    :param deg: Degree of K3
    :param n: Length of Hilbert scheme
    :param list_of_addtional: List of (mVec, mVec2, str: color).
     For each mVec draws additionally walls v(mVec) = v(mVec2) in color.
    """
    name = os.path.join(STORAGE_PATH_DIR, f"Degree{deg}Hilb{n}.walls")
    if not os.path.exists(name):
        print(f"Degree {deg} Hilb[{n}] not yet computed.")
        return
    file = open(name, "r")
    Hilb, walls_list = read_walls_file(file)
    file.close()

    k3 = Hilb.get_k3()
    muk = Hilb.get_mukai()
    print(str(Hilb) + " of H^2 = " + str(k3.get_deg()) + "\n")

    list_of_circles = []
    for w in walls_list:
        list_of_circles.append((w.get_vec(), muk, "black"))
    if list_of_addtional is not None:
        list_of_circles.extend(list_of_addtional)

    ModuliSpace.plot_graph(list_of_circles)


def show_all(list_of_deg_n: list[tuple[int, int]], rational: bool = False, show_walls: bool = True) -> None:
    """
    For a list of tuples (degree, n) show all walls for Hilb^n(K3 of degree)

    :param list_of_deg_n: List of tuples (degree, n)
    :param rational: Show the walls as fractions or decimal
    :param show_walls: Show walls (True) or only movable divisors (False)
    """
    for i in list_of_deg_n:
        deg, n = i
        name = os.path.join(STORAGE_PATH_DIR, f"Degree{deg}Hilb{n}.walls")
        if not os.path.exists(name):
            print(f"Degree {deg} Hilb[{n}] not yet computed.")
            continue
        file = open(name, "r")
        hilb, walls_list = read_walls_file(file)
        file.close()

        k3 = hilb.get_k3()
        print(str(hilb) + " of H^2 = " + str(k3.get_deg()) + "\n")

        basis_1 = mVec(k3, Rational(-1, 1), Rational(0, 1), Rational(-n+1, 1))
        basis_2 = mVec(k3, Rational(0, 1), Rational(-1, 1), Rational(0, 1))
        hilb.print_movable_divisors(walls_list, basis_1, basis_2, rational=rational, name_basis_1="B", name_basis_2="H")
        print(" ")
        if not show_walls:
            continue
        for w in walls_list:
            print(w)


HELP_STRING = """-h Help\n -s [a,b] """


def main_interactive() -> None:
    """
    The function that carries the interactive interface.
    """
    rational = True
    show_walls = True
    print("Initiated Hilbert scheme program with rational coefficients.")
    while True:
        print("""Type show|compute|showcompute|plot|rational|showwalls|depth|help|quit (s|c|sc|p|r|sw|d|h|q)""")
        inp = input()
        if inp in ["r", "ra", "rat", "rational"]:
            rational = not rational
            print(f"Rational toggled to {rational}")
            continue
        elif inp in ["w", "sw", "walls", "showwalls"]:
            show_walls = not show_walls
            print(f"Show_walls toggled to {show_walls}")
            continue
        elif inp in ["d", "depth"]:
            string = f"""Change search depth from \
{modulispace.modulispace.MAX_DEPTH_DIVISORIAL},\
{modulispace.modulispace.MAX_DEPTH_SPHERICAL_FLOPS},\
{modulispace.modulispace.MAX_DEPTH_SUM_FLOPS}\
 to user specified. Be aware that depth > 0 and depth = 1\
  can be risky and may lead to errors even in low degrees. It is advised to keep divisorial depth >= 2."""
            print(string)
            modulispace.modulispace.MAX_DEPTH_DIVISORIAL = int(input("Depth divisorial:    "))
            modulispace.modulispace.MAX_DEPTH_SPHERICAL_FLOPS = int(input("Depth sperical flop: "))
            modulispace.modulispace.MAX_DEPTH_SUM_FLOPS = int(input("Depth sum flop:      "))
            continue
        elif inp in ["", "q", "qu", "qui", "quit"]:
            print("Quitting...")
            break
        elif inp in ["h", "help"]:
            print(HELP_STRING_INTERACTIVE)
            continue

        b = inp in ["p", "pl", "plo", "plot"]
        list_d = get_degrees(only_one=b)
        if list_d is None:
            continue
        list_n = get_numbers(only_one=b)
        if list_n is None:
            continue
        list_all = [(d, n) for d in list_d for n in list_n]
        if inp in ["s", "sh", "sho", "show"]:
            show_all(list_all, rational, show_walls=show_walls)
        elif inp in ["c", "compute"]:
            for (d, n) in list_all:
                hilb, walls_list = compute_hilb(d, n)
                safe_hilb(hilb, walls_list, force_overwrite=True)
        elif inp in ["sc", "showcom", "showcompute"]:
            for (d, n) in list_all:
                if not exists_hilb(d, n):
                    hilb, walls_list = compute_hilb(d, n)
                    safe_hilb(hilb, walls_list)
                show_all([(d, n)], rational, show_walls=show_walls)
        elif inp in ["p", "pl", "plo", "plot"]:
            (d, n) = list_all[0]
            inp_add = None
            k3 = K3Surface(d)
            add_list = []
            while inp_add != "":
                inp_add = input("Additional vector: rk div xi (color): ")
                if inp_add == "":
                    break
                inp_add = inp_add.split()
                muk = mVec(k3, Rational(1, 1), Rational(0, 1), Rational(-n+1, 1))
                vec = mVec(k3, Rational(int(inp_add[0])), Rational(int(inp_add[1])), Rational(int(inp_add[2])))
                color = "red"
                if len(inp_add) > 3:
                    color = inp_add[3]
                add_list.append((muk, vec, color))
            show_plot(d, n, add_list)
        else:
            print(f"Keyword '{inp}'not found")
    return


def get_numbers(only_one: bool = False) -> Optional[list[int]]:
    """
    Asks user for numbers. Returns all numbers that the user wrote. Allows for rows, i.e. 2 - 4 is {2,3,4}.
    Also ".","..", "..." is allowed to inidcate a row.

    :param only_one:  Only one number is allowed as response
    :return: List of intended numbers.
    """
    if only_one:
        print("Number of points: ")
    else:
        print("Numbers of points: ")
    inp = input()
    inp = inp.split(" ")
    if len(inp) > 1 and only_one:
        print("Only one number allowed")
        return None
    list_numbers = []

    d_bool = False
    for s in inp:
        if s == "":
            continue
        if s in [".", "..", "...", "-"]:
            d_bool = True
            continue
        integer = int(s)
        if not d_bool:
            list_numbers.append(integer)
        else:
            last = list_numbers[-1]
            list_numbers.extend(list(range(last + 1, integer + 1)))
            d_bool = False
    return list_numbers


def get_degrees(only_one: bool = False) -> Optional[list[int]]:
    """
        Asks user for numbers. Returns all even(!) numbers that the user wrote. Allows for rows, i.e. 2 - 4 is {2,3,4}.
        Also ".","..", "..." is allowed to inidcate a row.
        Returns None if only_one is violated.

        :param only_one:  Only one number is allowed as response
        :return: List of intended numbers.
        """
    if only_one:
        print("Degree of K3: ")
    else:
        print("Degrees of K3s: ")
    inp = input()
    inp = inp.split(" ")
    if len(inp) > 1 and only_one:
        print("Only one degree allowed")
        return None
    list_degrees = []

    d_bool = False # if user inputted a row.
    for s in inp:
        if s == "":
            continue
        if s in [".", "..", "...", "-"]:
            d_bool = True
            continue
        integer = int(s)
        if not d_bool:
            list_degrees.append(integer)
        else:
            last = list_degrees[-1]
            list_degrees.extend(list(range(last + 2, integer + 1, 2)))
            d_bool = False
    return list_degrees


def thres_wall(wall: Wall) -> float:
    """
    Returns div/r if div>0 and -div/r otherwise.

    :param wall: Wall
    """
    divisor = wall.get_movable_divisor()
    r = divisor.get_rank().evalf()
    d = divisor.get_divisor().evalf()
    if d == 0:
        if r < 0:
            return -modulispace.modulispace._INFINITY
        return modulispace.modulispace._INFINITY
    if d < 0:
        r = -d/r
    else:
        r = r/d

    return r


if __name__ == "__main__":

    main_interactive()

