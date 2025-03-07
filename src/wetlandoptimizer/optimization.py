import numpy as np # type: ignore
import cma  # type: ignore
from wetlandoptimizer import treatment

class Optimizer_French_VF:
    """
    A class to perform CMA-ES optimization for the French vertical flow wetland.

    Attributes
    ----------
    treatment_train : Treatment_Train
        The treatment train to optimize.

    Methods
    -------
    Optimize():
        Performs the optimization.
    """
    def __init__(self, treatment_train):
        """
        Constructs the optimizer for French vertical flow wetland object.

        Parameters
        ----------
        treatment_train : Treatment_Train
            The treatment train to optimize.
        """
        self.treatment_train = treatment_train

    def Optimize(self):
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
        
        for process in self.treatment_train.pathway:
            bounds_min.extend([float(process.Xmin), float(process.Zmin), float(process.Amin)])
            bounds_max.extend([float(process.Xmax), float(process.Zmax), float(process.Amax)])

        bounds = [bounds_min, bounds_max]

        x0 = [(bounds[0][i] + bounds[1][i]) / 2 for i in range(len(bounds[0]))]

        sigma0 = 1e-2
        opts = {'seed': random_seed, 'bounds': bounds, 'popsize': 100}

        es = cma.CMAEvolutionStrategy(x0, sigma0, opts)

        while not es.stop():
            candidate_solutions = es.ask()
            fitness_values = [self.treatment_train.Single_Objective_Function(s, self.treatment_train.Cin, self.treatment_train.Cobj, self.treatment_train.Q) for s in candidate_solutions]
            es.tell(candidate_solutions, fitness_values)

        best_solution_cma = es.result.xbest
        return best_solution_cma

######################################################################################

class Optimizer_Gobal_Generation:
    """
    A class to perform CMA-ES optimization for the French vertical flow wetland.

    Attributes
    ----------
    treatment_train : Treatment_Train
        The treatment train to optimize.

    Methods
    -------
    Optimize():
        Performs the optimization for a single treatment train.
    Optimize_Best_Pathway():
        Generate the possible combinations and call the Optimize() function for each generated treatment train.
    """
    def __init__(self, treatment_train, stages_max, files_max):
        """
        Constructs the optimizer for a treatment train object.

        Parameters
        ----------
        treatment_train : Treatment_Train
            The treatment train to optimize.
        stages_max : float
            Maximum value for the number of stages in series.
        files_max : float
            Maximum value for the number of files in parallel.
        """
        self.treatment_train = treatment_train
        self.stages_max = stages_max
        self.files_max = files_max

    def Optimize(self, pathway):
        """
        Performs the optimization.

        Parameters
        ----------
        pathway : Pathway
            The pathway of the treatment train to optimize.        

        Returns
        -------
        best_solution_cma : list
            Best solution found by CMA-ES.
        """
        random_seed = np.random.seed()
        # total_decision_vars = sum(3 for process in pathway)
        
        bounds_min = []
        bounds_max = []
        
        # for process in self.treatment_train.pathway:
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
        
        self.treatment_train.pathway = pathway
        
        while not es.stop():
            candidate_solutions = es.ask()
            # print(self.treatment_train.pathway)
            # print(self.treatment_train.pathway)
            fitness_values = [self.treatment_train.Single_Objective_Function(s, self.treatment_train.Cin, self.treatment_train.Cobj, self.treatment_train.Q) for s in candidate_solutions]
            es.tell(candidate_solutions, fitness_values)

        best_solution_cma = es.result.xbest
        # print("test2",best_solution_cma)

        # es = None
        return best_solution_cma
      
    def Optimize_Best_Pathway(self):
        """
        Generate the possible combinations and call the Optimize() function for each generated treatment train.

        Returns
        -------
        sorted_results : list
            List of optimized and possible combinations.
        """
        pathway = treatment.Pathway(self.stages_max, self.files_max)
        pathways_possibles = pathway.Possible_Combinations()
        print(pathways_possibles)

        best_results = []
        rejected_results = []

        for pathway_combination in pathways_possibles:
            solution = self.Optimize(pathway_combination)

            if solution is None:
                print("No suitable solutions for this combination of Cin / Cobj")
                rejected_results.append([pathway_combination, None])
                continue

            self.treatment_train.pathway = pathway_combination
            volume = self.treatment_train.Total_Volume_Function(solution, self.treatment_train.Q)

            constraints_TSSout = self.treatment_train.Create_Constraints_TSSout(solution, self.treatment_train.Cobj, self.treatment_train.Q)
            constraints_BODout = self.treatment_train.Create_Constraints_BODout(solution, self.treatment_train.Cobj, self.treatment_train.Q)
            constraints_TKNout = self.treatment_train.Create_Constraints_TKNout(solution, self.treatment_train.Cobj, self.treatment_train.Q)
            constraints_CODout = self.treatment_train.Create_Constraints_CODout(solution, self.treatment_train.Cobj, self.treatment_train.Q)
            constraints_NO3out = self.treatment_train.Create_Constraints_NO3out(solution, self.treatment_train.Cobj, self.treatment_train.Q)
            constraints_TNout = self.treatment_train.Create_Constraints_TNout(solution, self.treatment_train.Cobj, self.treatment_train.Q)            

            constraints_met = all(c >= 0 for c in constraints_TSSout) and \
                            all(c >= 0 for c in constraints_BODout) and \
                            all(c >= 0 for c in constraints_TKNout) and \
                            all(c >= 0 for c in constraints_CODout) and \
                            all(c >= 0 for c in constraints_NO3out) and \
                            all(c >= 0 for c in constraints_TNout)

            if constraints_met:
                best_results.append([pathway_combination, volume, solution])
            if not constraints_met:
                rejected_results.append([pathway_combination, solution])

        rejected_names = [
            [type(process).__name__ for process in pathway_combination]
            for pathway_combination, _ in rejected_results
        ]
        print(f"Rejected solutions: {len(rejected_results)} ({', '.join(['-'.join(names) for names in rejected_names])})")
        print("")

        sorted_results = sorted(best_results, key=lambda x: x[1])

        return sorted_results
