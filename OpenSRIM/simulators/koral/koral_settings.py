class KORALSettings():
    def __init__(self,
                 nr_iterations: int = 1,
                 integration_method: str = 'LSODA',
                 rtol: float = 1e-6,
                 atol: float = 1e-6):
        """Creates an object containing all settings
        for KORAL calculation

        Args:
            nr_iterations (int, optional): Number of iterations of KORAL.
                Defaults to 1.
            integration_method (str, optional): Method of integration.
                Defaults to 'RK45'.
            rtol (float, optional): Total relative integration error.
                Defaults to 1e-6.
            atol (float, optional): Total absolute integration error.
                Defaults to 1e-6.
        """

        self.nr_interations = nr_iterations
        self.integration_method = integration_method
        self.rtol = rtol
        self.atol = atol

        # TODO: Werte validieren