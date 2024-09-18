import numpy as np # type: ignore
import cma  # type: ignore

class CMAES:
    """
    A class to perform CMA-ES optimization on the treatment chain.

    Attributes
    ----------
    treatment_chain : Treatment_Chain
        The treatment chain to optimize.

    Methods
    -------
    optimize():
        Performs the optimization.
    """
    def __init__(self, treatment_chain):
        """
        Constructs the CMAES optimizer object.

        Parameters
        ----------
        treatment_chain : Treatment_Chain
            The treatment chain to optimize.
        """
        self.treatment_chain = treatment_chain
    
    def optimize(self):
        """
        Performs the optimization.

        Returns
        -------
        best_solution_cma : list
            Best solution found by CMA-ES.
        """
        random_seed = np.random.seed()
        
        bounds_min = []
        bounds_max = []
        
        for process in self.treatment_chain.pathway:
            bounds_min.extend([process.Xmin, process.Zmin])
            bounds_max.extend([process.Xmax, process.Zmax])
          
        bounds = [bounds_min, bounds_max]

        x0 = [(bounds[0][i] + bounds[1][i]) / 2 for i in range(len(bounds[0]))]
        
        sigma0 = 1e-3
        opts = {'seed': random_seed, 'bounds': bounds, 'popsize': 100}

        es = cma.CMAEvolutionStrategy(x0, sigma0, opts)

        while not es.stop():
            candidate_solutions = es.ask()
            fitness_values = [self.treatment_chain.Single_Objective_Function(s, self.treatment_chain.Cin, self.treatment_chain.Cobj, self.treatment_chain.Q) for s in candidate_solutions]
            es.tell(candidate_solutions, fitness_values)

        best_solution_cma = es.result.xbest

        return best_solution_cma
