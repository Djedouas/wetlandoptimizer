import numpy as np # type: ignore
import threading
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
      
    def process_combination(self, pathway_combination, best_results_lock, rejected_results_lock, best_results, rejected_results):
        solution = self.Optimize(pathway_combination)

        if solution is None:
            with rejected_results_lock:
                print(f"No suitable solutions for this combination of Cin / Cobj: {pathway_combination}")
                rejected_results.append([pathway_combination, None])
            return

        self.treatment_train.pathway = pathway_combination
        volume = self.treatment_train.Total_Volume_Function(solution, self.treatment_train.Q)

        # Créer les différentes contraintes
        constraints_list = [
            self.treatment_train.Create_Constraints_TSSout(solution, self.treatment_train.Cobj, self.treatment_train.Q),
            self.treatment_train.Create_Constraints_BODout(solution, self.treatment_train.Cobj, self.treatment_train.Q),
            self.treatment_train.Create_Constraints_TKNout(solution, self.treatment_train.Cobj, self.treatment_train.Q),
            self.treatment_train.Create_Constraints_CODout(solution, self.treatment_train.Cobj, self.treatment_train.Q),
            self.treatment_train.Create_Constraints_NO3out(solution, self.treatment_train.Cobj, self.treatment_train.Q),
            self.treatment_train.Create_Constraints_TNout(solution, self.treatment_train.Cobj, self.treatment_train.Q)
        ]

        # Vérification des contraintes
        constraints_met = all(all(c >= 0 for c in constraint) for constraint in constraints_list)

        if constraints_met:
            with best_results_lock:
                best_results.append([pathway_combination, volume, solution])
        else:
            with rejected_results_lock:
                rejected_results.append([pathway_combination, solution])
                # Optionnel : Ajouter des informations sur quelles contraintes ont échoué
                for i, constraint in enumerate(constraints_list):
                    if any(c < 0 for c in constraint):
                        print(f"Constraint {i+1} violated for pathway {pathway_combination}")

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

        threads = []
        best_results = []
        rejected_results = []
        best_results_lock = threading.Lock()
        rejected_results_lock = threading.Lock()

        for pathway_combination in pathways_possibles:
            # Créer un thread pour chaque combinaison
            thread = threading.Thread(target=self.process_combination,
                                      args=(pathway_combination, best_results_lock, rejected_results_lock, best_results, rejected_results)
                                      )
            threads.append(thread)
            thread.start()

        # Attendre que tous les threads aient terminé
        for thread in threads:
            thread.join()

        rejected_names = [
            [type(process).__name__ for process in pathway_combination]
            for pathway_combination, _ in rejected_results
        ]
        print(f"Rejected solutions: {len(rejected_results)} ({', '.join(['-'.join(names) for names in rejected_names])})")
        print("")

        sorted_results = sorted(best_results, key=lambda x: x[1])

        return sorted_results
