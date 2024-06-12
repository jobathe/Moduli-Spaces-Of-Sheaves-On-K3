from modulispace import ModuliSpace
from modulispace import mVec
from modulispace import K3Surface
from modulispace import write_walls_file, read_walls_file


def compute_hilbert_scheme():
    """In this function we compute Hilb^3 of a K3 of degree 2"""

    print("Compute walls of Hilb^[3] of a K3 of degree 2.")

    # Degree 2 K3 surface
    k3 = K3Surface(2)

    # Hilbert scheme of 3 points = M(1,0,-2)
    hilbert = ModuliSpace(k3, mVec(k3, 1, 0, -2))

    # Specify ample vector 100H^[3] - B
    ample = mVec(k3, 1, -100, 2)

    # Compute walls
    walls = hilbert.compute_walls(ample)

    # The divisor B = E/2, where E is the divisor of non-red. schemes
    basis_1 = mVec(k3, -1, 0, -2)

    # The divisor H^[3]
    basis_2 = mVec(k3, 0, -1, 0)

    # print movable divisors and the contracted effective divisors
    hilbert.print_movable_divisors(walls, basis_1, basis_2, name_basis_1="B", name_basis_2="H", print_contracted=True)

    # print the walls
    for w in walls:
        print(w)

    with open("examples.walls", "w") as file:
        # write the walls to the file "examples.walls"
        write_walls_file(file, hilbert, walls, ample=ample)

    return


def open_walls_file():
    """This function opens the file 'examples.walls' and prints the Hilbert scheme and the walls."""
    with open("examples.walls", "r") as file:
        hilbert, walls = read_walls_file(file)

    print("We just loaded the file 'examples.wall'. It consists of the following: ")

    print(hilbert)
    print(hilbert.get_k3())

    for w in walls:
        print(w)


compute_hilbert_scheme()
open_walls_file()
