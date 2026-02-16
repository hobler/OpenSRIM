class KORALInput():
    def __init__(self,
                 method: str,
                 z_ion: int,
                 m_ion: float,
                 z_target: list[int],
                 m_target: list[float],
                 d_target: float,
                 s_e_f: list[float],
                 start_energy: int,
                 stop_energy: int,
                 nr_values: int):    
        """Creates an object containing all input paramerters
        for KORAL calculation

        Args:
            method (str): Calculation method.
                Supported values: 'ZBL', 'NLH'
            z_ion (int): Atomic number of ion
            m_ion (float): Relative atomic mass of ion
            z_target (list[int]): Atomic numbers of targets
            m_target (list[float]): Relative atomic masses of targets
            d_target (float): Element density of target (1 / A^3)
            s_e_f (list[float]): Stopping power weight factors
            start_energy (int): Calculation starts at this energy (eV)
            stop_energy (int): Calculation stops at this energy (eV)
            nr_values (int): Number of calculated result values between
                start_energy and stop_energy
        """
        self.method = method
        self.z_ion = z_ion
        self.m_ion = m_ion
        self.z_target = z_target
        self.m_target = m_target
        self.d_target = d_target
        self.s_e_f = s_e_f
        self.start_energy = start_energy
        self.stop_energy = stop_energy
        self.nr_values = nr_values
        # TODO: Woher kommen die S_e Werte? Müssen die aus den Files geladen werden?
        # TODO: Selbes gilt für die NLH Parameter für S_n, Q_n und V