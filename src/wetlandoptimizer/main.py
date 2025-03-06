import treatment
import optimization
import yaml
import sys
import os
import numpy as np

def COD_Fractionation(pollutants_in):
    """
    Fractionates the total COD (CODt) in three parts: COD dissolved biodegradable part (CODdb), COD dissolved inert part (CODdi) and COD particulate part (CODp).

    Parameters
    ----------
    polluants_in : list
        Unfractionnated input concentrations ([TSS]in : polluants_in[0] (gTSS/m3), [BOD5]in : polluants_in[1] (gO2/m3), [TKN]in : polluants_in[2] (gTKN/m3), [CODt]in : polluants_in[3] (gO2/m2)).

    Returns
    -------
    fractionated_polluants_in : list
        Fractionated input concentrations ([TSS]in : fractionnated_polluants_in[0] (gTSS/m3), [BOD5]in : fractionnated_polluants_in[1] (gO2/m3), [TKN]in : fractionnated_polluants_in[2] (gTKN/m3), [CODdb]in : fractionnated_polluants_in[3] (gO2/m2), [CODdi]in : fractionnated_polluants_in[4] (gO2/m2), [CODp]in : fractionnated_polluants_in[5] (gO2/m2)).
    """  
    if pollutants_in[0] > 0.75 * pollutants_in[3] :
        print("erreur")
        fractionated_pollutants_in = [0,0,0,0,0]
    else :  
        P = 1.1 * pollutants_in[0]
        if 0.04 * pollutants_in[3] >= 30 :
            Si = 30.0
        else :
            Si = 0.04 * pollutants_in[3]
        Sb = pollutants_in[3] - P - Si
        fractionated_pollutants_in = [pollutants_in[0],pollutants_in[1],pollutants_in[2],Sb,Si,P,pollutants_in[4]]
    return fractionated_pollutants_in

################################################################################

def Results_French_VF(Cin, Cobj, Q, climate) :
    """
    Optimizes the French vertical flow wetland and prints the results for the user.

    Parameters
    ----------
    Cin : list
        Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
    Cobj : list
        Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
    Q : float
        Flow rate (m3/day).
    climate : str
        Climate type, either "temperate" or "tropical", determining environmental conditions affecting the treatment process.
    """
    base_directory = r'C:\Users\zoe.legeai\Documents\Source\wetlandoptimizer\src\wetlandoptimizer'
    sys.path.append(base_directory)

    config_path = os.path.join(base_directory, 'config.yaml')
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"No such file or directory: '{config_path}'")
    
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    if climate == "Tropical":
        config['VdNS1']['Nb_parallel'] = 2

    vdns1 = treatment.VdNS1(**config['VdNS1'])
    vdns2 = treatment.VdNS2(**config['VdNS2'])
    pathway_french_TW = [vdns1, vdns2]

    print("Optimizing with CMA-ES...")
    print("")

    treatment_train_french_TW = treatment.Treatment_Train(pathway_french_TW, [], Cin, Cobj, Q)
    optimizer = optimization.Optimizer_French_VF(treatment_train_french_TW)
    best_solution_cma = optimizer.Optimize()

    treatment_train_french_TW_opti = treatment.Treatment_Train(pathway_french_TW, best_solution_cma, Cin, Cobj, Q)
    output_function_values = treatment_train_french_TW_opti.Output_Function(best_solution_cma, Q)
    
    print("")
    print("Best solution:", [round(value, 2) for value in best_solution_cma])
    print("")
    print("Sizing:")
    print("")
    print("Total surface area 1st floor :",round(Q/best_solution_cma[0]*pathway_french_TW[0].Nb_parallel,2),"m2")
    print("Depth 1st floor :",round(best_solution_cma[1],2),"m")
    print("---")
    if len(pathway_french_TW) == 2 :
        print("Total surface area 2nd floor :",round(Q/best_solution_cma[3]*pathway_french_TW[1].Nb_parallel,2),"m2")
        print("Depth 2nd floor :",round(best_solution_cma[4],2),"m")
    print("---")
    print("Total volume :",round(treatment_train_french_TW_opti.Total_Volume_Function(best_solution_cma, Q),2),"m3")
    print("Total surface area :",round(treatment_train_french_TW_opti.Total_Surface_Area_Function(best_solution_cma, Q),2),"m2")
    print("")
    print("Outlet concentration:")
    print("")
    print("TSS (mgTSS/L):", round(output_function_values[0], 2))
    print("BOD5 (mgO2/L):", round(output_function_values[1], 2))
    print("TKN (mgTKN/L):", round(output_function_values[2], 2))
    print("COD (mgO2/L):", round(output_function_values[3] + output_function_values[4] + output_function_values[5], 2))
    print("NO3 (mgNO3/L):", round(output_function_values[6], 2))
    print("TN (mgNO3/L):", round(output_function_values[2]+output_function_values[6], 2))

    print("")
    print("Checking constraints values for best solution found by CMA-ES...")
    print("---")
    output = Cin
    for index, process in enumerate(pathway_french_TW):
        V_values = best_solution_cma[index * 3: (index + 1) * 3]     
        constraint_TSS_pc = 100 * output[0] * V_values[0] / process.Lim_TSS 
        print("TSS loading contraint stage", index+1, ":", round(constraint_TSS_pc,2), "%")
        constraint_BOD_pc = 100 * output[1] * V_values[0] / process.Lim_BOD
        print("BOD loading contraint stage", index+1, ":", round(constraint_BOD_pc,2), "%")
        constraint_TKN_pc = 100 * output[2] * V_values[0] / process.Lim_TKN
        print("TKN loading contraint stage", index+1, ":", round(constraint_TKN_pc,2), "%")
        constraint_COD_pc = 100 * (output[3]+output[4]+output[5]) * V_values[0] / process.Lim_COD
        print("COD loading contraint stage", index+1, ":", round(constraint_COD_pc,2), "%")
        print("")
        constraint_hydraulic_pc = 100 * (V_values[0] - process.Xmin) / (process.Xmax - process.Xmin)
        print("Hydraulic loading constraint stage", index+1, ":", round(constraint_hydraulic_pc,2), "%")
        print("---")
        output = process.Reduction_Function(V_values, output, Q)

    print("Outlet TSS deviation:",round(-(Cobj[0] - output_function_values[0]),2),"mgTSS/L")
    print("Outlet BOD5 deviation:",round(-(Cobj[1] - output_function_values[1]),2),"mgO2/L")
    if Cobj[2] != None:        
        print("Outlet TKN deviation:",round(-(Cobj[2] - output_function_values[2]),2),"mgTKN/L")
    print("Outlet COD deviation:",round(-(Cobj[3] - (output_function_values[3]+output_function_values[4]+output_function_values[5])),2),"mgO2/L")
    if Cobj[4] != None :
        print("Outlet NO3 deviation:",round(-(Cobj[4] - output_function_values[6]),2),"mgNO3/L")
    if Cobj[5] != None:    
        print("Outlet TN deviation:",round(-(Cobj[5] - (output_function_values[2]+output_function_values[6])),2),"mgO2/L")
      
################################################################################

def get_ordinal_suffix(number):
    """Returns the correct ordinal suffix for a number."""
    if 10 <= number % 100 <= 20:  # Gestion des exceptions comme 11th, 12th, 13th
        return "th"
    elif number % 10 == 1:
        return "st"
    elif number % 10 == 2:
        return "nd"
    elif number % 10 == 3:
        return "rd"
    else:
        return "th"
    
################################################################################

def Results_Global_Generation(Cin, Cobj, Q, stages_max, files_max, climate):
    """
    Optimizes and ranks the treatment trains generated exploring all the process combinations, ranks them by volume and prints the results for the user.

    Parameters
    ----------
    Cin : list
        Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
    Cobj : list
        Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
    Q : float
        Flow rate (m3/day).
    climate : str
        Climate type, either "temperate" or "tropical", determining environmental conditions affecting the treatment process.
    stages_max : float
        Maximum value for the number of stages in series.
    files_max : float
        Maximum value for the number of files in parallel.
    """
    base_directory = r'C:\Users\zoe.legeai\Documents\Source\caribsan-model\CARIBSANcopy\src\CARIBSANcopy'
    sys.path.append(base_directory)

    config_path = os.path.join(base_directory, 'config.yaml')
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"No such file or directory: '{config_path}'")
    
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    if climate == "Tropical":
        config['VdNS1']['Nb_parallel'] = 2
        config['VdNSS']['Nb_parallel'] = 2

    process_mapping = {}
    subclasses = treatment.Pathway.Get_Subclasses(treatment.Process)
    for subclass in subclasses:
        process_mapping[subclass.__name__] = subclass

    print("Optimizing with CMA-ES...")
    print("")

    all_results = []
    for current_stages_max in range(1, stages_max + 1):
        print(f"Optimizing for stages_max = {current_stages_max}...")
        treatment_train_MS_TW = treatment.Treatment_Train([], [], Cin, Cobj, Q)
        optimizer = optimization.Optimizer_Gobal_Generation(treatment_train_MS_TW, current_stages_max, files_max)
        sorted_results = optimizer.Optimize_Best_Pathway()

        for result in sorted_results:
            pathway_cleaned = [type(obj).__name__ for obj in result[0]]
            solution = result[2]
            total_volume = result[1]

            if not any(
                existing_result['pathway'] == pathway_cleaned and 
                np.allclose(existing_result['solution'], solution, atol=1e-3)
                for existing_result in all_results):
                all_results.append({
                    'stages_max': current_stages_max,
                    'pathway': pathway_cleaned,
                    'solution': solution,
                    'volume': total_volume
                })

    all_results = sorted(all_results, key=lambda x: x['volume'])
    if not all_results:
        print("No suitable solutions for this combination of Cin / Cobj")

    for position, result in enumerate(all_results):
        print("")
        print("################################################################################")
        print("")
        print(f"stages_max: {result['stages_max']}, Pathway at position: {(position + 1)}")
        print("Pathway:", result['pathway'])
        print("")
        print("Solution:", [round(value, 2) for value in result['solution']])
        print("")

        process_instances = [
            process_mapping[process_name](**config[process_name])
            for process_name in result['pathway']
        ]

        treatment_train_MS_TW_opti = treatment.Treatment_Train(
            process_instances,
            result['solution'],
            Cin, Cobj, Q
        )

        output_function_values = treatment_train_MS_TW_opti.Output_Function(result['solution'], Q)
        total_volume = treatment_train_MS_TW_opti.Total_Volume_Function(result['solution'], Q)

        print("Sizing:")
        print("")

        for stage_index, pathway in enumerate(result['pathway']):
            process_class = process_mapping[pathway]
            process_instance = process_class(**config[pathway]) 

            surface_area = round(Q / result['solution'][stage_index * 3] * process_instance.Nb_parallel, 2)
            depth_unsat = round(result['solution'][stage_index * 3 + 1], 2)
            depth_sat = round(result['solution'][stage_index * 3 + 2], 2)

            ordinal_suffix = get_ordinal_suffix(stage_index + 1)

            print(f"Total surface area {stage_index + 1}{ordinal_suffix} floor :", surface_area, "m²")
            print(f"Unsaturated depth {stage_index + 1}{ordinal_suffix} floor :", depth_unsat, "m")
            print(f"Saturated depth {stage_index + 1}{ordinal_suffix} floor :", depth_sat, "m")
            print("---")

        print("Total volume:", round(treatment_train_MS_TW_opti.Total_Volume_Function(result['solution'], Q), 2), "m³")
        print("Total surface area:", round(treatment_train_MS_TW_opti.Total_Surface_Area_Function(result['solution'], Q), 2), "m²")

        print("")
        print("Outlet concentration:")
        print("")
        print("TSS:", round(output_function_values[0], 2), "mgTSS/L")
        print("BOD5:", round(output_function_values[1], 2), "mgO2/L")
        print("TKN:", round(output_function_values[2], 2), "mgTKN/L")
        print("COD:", round(output_function_values[3] + output_function_values[4] + output_function_values[5], 2), "mgO2/L")
        print("NO3:", round(output_function_values[6], 2), "mgNO3/L")
        print("TN:", round(output_function_values[2]+output_function_values[6], 2), "mgN/L")

        print("")
        print("Checking constraints values for solution found by CMA-ES...")
        print("")
        output = Cin

        for index, process_name in enumerate(result['pathway']):
            process_class = process_mapping[process_name]
            process_instance = process_class(**config[process_name])
            V_values = result['solution'][index * 3: (index + 1) * 3]

            constraint_TSS_pc = 100 * output[0] * V_values[0] / process_instance.Lim_TSS
            print(f"TSS loading contraint stage {index + 1}: {round(constraint_TSS_pc, 2)} %")
            constraint_BOD_pc = 100 * output[1] * V_values[0] / process_instance.Lim_BOD
            print(f"BOD loading contraint stage {index + 1}: {round(constraint_BOD_pc, 2)} %")
            constraint_TKN_pc = 100 * output[2] * V_values[0] / process_instance.Lim_TKN
            print(f"TKN loading contraint stage {index + 1}: {round(constraint_TKN_pc, 2)} %")
            constraint_COD_pc = 100 * (output[3] + output[4] + output[5]) * V_values[0] / process_instance.Lim_COD
            print(f"COD loading contraint stage {index + 1}: {round(constraint_COD_pc, 2)} %")
            print("---")
            
            output = process_instance.Reduction_Function(V_values, output, Q)

        print("Outlet TSS deviation:", round(-(Cobj[0] - output_function_values[0]), 2), "mgTSS/L")
        print("Outlet BOD5 deviation:", round(-(Cobj[1] - output_function_values[1]), 2), "mgO2/L")
        if Cobj[2] != None:        
            print("Outlet TKN deviation:", round(-(Cobj[2] - output_function_values[2]), 2), "mgTKN/L")
        print("Outlet COD deviation:", round(-(Cobj[3] - (output_function_values[3] + output_function_values[4] + output_function_values[5])), 2), "mgO2/L")
        if Cobj[4] != None:
            print("Outlet NO3 deviation:", round(-(Cobj[4] - output_function_values[6]), 2), "mgNO3/L")
        if Cobj[5] != None:
            print("Outlet TN deviation:", round(-(Cobj[5] - (output_function_values[2] + output_function_values[6])), 2), "mgN/L")
