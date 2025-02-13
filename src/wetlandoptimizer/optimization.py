import numpy as np # type: ignore
import cma  # type: ignore
import treatment

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
            bounds_min.extend([float(process.Xmin), float(process.Zmin), float(process.Amin)])
            bounds_max.extend([float(process.Xmax), float(process.Zmax), float(process.Amax)])

        bounds = [bounds_min, bounds_max]

        x0 = [(bounds[0][i] + bounds[1][i]) / 2 for i in range(len(bounds[0]))]

        sigma0 = 1e-2
        opts = {'seed': random_seed, 'bounds': bounds, 'popsize': 100}

        es = cma.CMAEvolutionStrategy(x0, sigma0, opts)

        while not es.stop():
            candidate_solutions = es.ask()
            fitness_values = [self.treatment_chain.Single_Objective_Function(s, self.treatment_chain.Cin, self.treatment_chain.Cobj, self.treatment_chain.Q) for s in candidate_solutions]
            es.tell(candidate_solutions, fitness_values)

        best_solution_cma = es.result.xbest
        return best_solution_cma

######################################################################################

class Optimiseur5:

    # base_directory = r'C:\Users\zoe.legeai\Documents\Source\caribsan-model\CARIBSANcopy\src\CARIBSANcopy'
    # sys.path.append(base_directory)

    # config_path = os.path.join(base_directory, 'config.yaml')
    
    # if not os.path.exists(config_path):
    #     raise FileNotFoundError(f"No such file or directory: '{config_path}'")
    
    # with open(config_path, 'r') as file:
    #     config = yaml.safe_load(file)
    
    def __init__(self, treatment_chain, stages_max, files_max):
        self.treatment_chain = treatment_chain
        self.stages_max = stages_max
        self.files_max = files_max

    def Opti_CMAES5(self, pathway):
        random_seed = np.random.seed()
        # total_decision_vars = sum(3 for process in pathway)
        
        bounds_min = []
        bounds_max = []
        
        # for process in self.treatment_chain.pathway:
            # print(f"Process: {process}, Xmin: {process.Xmin}, Zmin: {process.Zmin}, Amin: {process.Amin}")

        for process in pathway:
            # print(f"Xmin: {process.Xmin}, Zmin: {process.Zmin}, Amin: {process.Amin}")
            bounds_min.extend([float(process.Xmin), float(process.Zmin), float(process.Amin)])
            bounds_max.extend([float(process.Xmax), float(process.Zmax), float(process.Amax)])

        bounds = [bounds_min, bounds_max]
        # print(bounds)

        x0 = [(bounds[0][i] + bounds[1][i]) / 2 for i in range(len(bounds[0]))]

        sigma0 = 1e-3
        opts = {'seed': random_seed, 'bounds': bounds, 'popsize': 100}

        es = cma.CMAEvolutionStrategy(x0, sigma0, opts)
        # print("test", es)
        
        self.treatment_chain.pathway = pathway
        
        while not es.stop():
            candidate_solutions = es.ask()
            # print(self.treatment_chain.pathway)
            # print(self.treatment_chain.pathway)
            fitness_values = [self.treatment_chain.Single_Objective_Function(s, self.treatment_chain.Cin, self.treatment_chain.Cobj, self.treatment_chain.Q) for s in candidate_solutions]
            es.tell(candidate_solutions, fitness_values)

        best_solution_cma = es.result.xbest
        # print("test2",best_solution_cma)

        # es = None
        return best_solution_cma
      

    def Opti_meilleur_pathway(self):
        pathway = treatment.Pathway(self.stages_max, self.files_max)
        pathways_possibles = pathway.Possible_Combinations()
        print(pathways_possibles)

        best_results = []
        rejected_results = []
        for pathway_combination in pathways_possibles:
            solution = self.Opti_CMAES5(pathway_combination)
            self.treatment_chain.pathway = pathway_combination
            volume = self.treatment_chain.Total_Volume_Function(solution, self.treatment_chain.Q)

            constraints_TSSout = self.treatment_chain.Create_Constraints_TSSout(solution, self.treatment_chain.Cobj, self.treatment_chain.Q)
            constraints_BODout = self.treatment_chain.Create_Constraints_BODout(solution, self.treatment_chain.Cobj, self.treatment_chain.Q)
            constraints_TKNout = self.treatment_chain.Create_Constraints_TKNout(solution, self.treatment_chain.Cobj, self.treatment_chain.Q)
            constraints_CODout = self.treatment_chain.Create_Constraints_CODout(solution, self.treatment_chain.Cobj, self.treatment_chain.Q)
            constraints_CODout = self.treatment_chain.Create_Constraints_NO3out(solution, self.treatment_chain.Cobj, self.treatment_chain.Q)

            constraints_met = all(c >= 0 for c in constraints_TSSout) and \
                            all(c >= 0 for c in constraints_BODout) and \
                            all(c >= 0 for c in constraints_TKNout) and \
                            all(c >= 0 for c in constraints_CODout)

            if constraints_met:
                best_results.append([pathway_combination, volume, solution])
            if not constraints_met:
                rejected_results.append([pathway_combination, solution])

        rejected_names = [
            [type(process).__name__ for process in pathway_combination]
            for pathway_combination, _ in rejected_results
        ]
        print("")
        print(f"Rejected solutions: {len(rejected_results)} ({', '.join(['-'.join(names) for names in rejected_names])})")

        sorted_results = sorted(best_results, key=lambda x: x[1])

        return sorted_results