from modulispace.lattice import Lattice


class SpecialClass:
    """Represents the vector that induces a Wall."""

    NON = -1  # Not a Wall
    BN = 0
    HC = 1
    LGU = 2
    FLSP = 3  # Flop sum
    FLSU = 4  # Flop spherical
    dict = {NON: "No Contraction", BN: "Brill-Noether", HC: "Hilbert-Chow",
            LGU: "Lie-Gieseker-Uhlenbeck", FLSP: "Spherical Flop", FLSU: "Sum Flop"}

    def __init__(self, moduli_space, vec):
        self.vec = vec
        self.moduli_space = moduli_space

    def get_vec(self):
        return self.vec

    def get_modulispace(self):
        return self.moduli_space

    def get_codim(self):
        return

    def get_type_as_str(self):
        return self.dict[self.get_type_num()]

    def get_type_num(self):
        """
        If vector induces a contraction returns the contraction type as SpecialClass.Contraction.
        Returns NON if the vector does not induce a contraction.
        """
        mukai = self.get_modulispace().get_mukai()
        vec = self.get_vec()

        lattice = Lattice(vec, mukai)
        if lattice.get_rank() == 1:
            return self.NON

        if not lattice.is_hyperbolic():
            return self.NON
        if vec*vec == 0:
            # isotropic vectors
            if vec*mukai == 1:
                return self.HC
            elif vec*mukai == 2:
                if mukai.get_rank() % 2 == 0 and mukai.get_divisor() % 2 == 0 and mukai.get_ch2() % 2 == 0:
                    return self.HC
                else:
                    return self.LGU
        elif vec*vec == -2:
            # first check flop as they are mutually exclusive to
            if vec * mukai > 0 and 2 * (vec * mukai) <= mukai * mukai:
                return self.FLSP
            # first check BN
            if lattice.contains_isotropic_element():
                # Check if it is a flop as it is not a divisorial contraction
                return self.NON
            if vec*mukai == 0:
                return self.BN

        if vec*vec >= 0 and vec*mukai > 0 and (mukai-vec)*(mukai-vec) >= 0 and (mukai-vec)*mukai > 0:
            return self.FLSU
        return self.NON
