class K3Surface:
    """
    Represents a K3 surface with Pic X = Z. The class only carries the information on the degree.
    """
    def __init__(self, d: int) -> None:
        if d % 2 != 0 or d <= 0:
            raise ValueError("Degree of a K3 surface must be an even integer > 0.")
        self.degree = d

    def get_deg(self) -> int:
        return self.degree

    def __str__(self) -> str:
        return f"K3 Surface with H.H = {self.get_deg()}"
