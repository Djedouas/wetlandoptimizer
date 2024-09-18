import treatment
import optimization
import yaml

def load_config(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def COD_Fractionation(polluants_in) :
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
    if polluants_in[0] > 0.75 * polluants_in[3] :
        
        print("erreur")
        fractionated_polluants_in = [0,0,0,0]
        
    else :  
        
        P = 1.1 * polluants_in[0]
    
        if 0.04 * polluants_in[3] >= 30 :
        
            Si = 30.0
        
        else :
        
            Si = 0.04 * polluants_in[3]
        
        Sb = polluants_in[3] - P - Si
    
        fractionated_polluants_in = [polluants_in[0],polluants_in[1],polluants_in[2],Sb,Si,P]
    
    return fractionated_polluants_in

################################################################################

def main(TSS_in,BOD5_in,TKN_in,COD_in,TSS_out,BOD5_out,TKN_out,COD_out,Q) :
    """
    Optimizes the treatment chain and prints the results for the user.

    Parameters
    ----------
    TSS_in : float
        Inlet TSS concentration (mgTSS/L).
    BOD5_in : float
        Inlet BOD5 concentration (mgO2/L).
    TKN_in : float
        Inlet TKN concentration (mgTKN/L).
    COD_in : float
        Inlet CODt concentration (mgO2/L).
    TSS_out : float
        Outlet TSS objective concentration (mgTSS/L).
    BOD5_out : float
        Outlet BOD5 objective concentration (mgO2/L).
    TKN_out : float
        Outlet TKN objective concentration (mgTKN/L).
    COD_out : float
        Outlet CODt objective concentration (mgO2/L).
    Q : float
        Flow rate (m3/day).
    """
    config_path = 'src/config.yaml'
    config = load_config(config_path)

    Cin = COD_Fractionation([TSS_in, BOD5_in, TKN_in, COD_in]) # g/m3
    Cobj = [TSS_out, BOD5_out, TKN_out, COD_out] # g/m3
  
    vdns1 = treatment.VdNS1(**config['VdNS1'])
    vdns2 = treatment.VdNS2(**config['VdNS2'])
    pathway_french_TW = [vdns1, vdns2]

    print("Optimizing with CMA-ES...")
    print("")
    treatment_chain_french_TW = treatment.Treatment_Chain(pathway_french_TW, [], Cin, Cobj, Q)
    optimiseur = optimization.CMAES(treatment_chain_french_TW)
    best_solution_cma = optimiseur.optimize()
    print("")
    print("Best solution:", [round(value, 2) for value in best_solution_cma])

    print("")
    print("Sizing:")
    print("")
    print("Total surface area 1st floor :",round(Q/best_solution_cma[0]*vdns1.Nb_parallel,2),"m2")
    print("Depth 1st floor :",round(best_solution_cma[1],2),"m")
    print("---")
    print("Total surface area 2nd floor :",round(Q/best_solution_cma[2]*vdns2.Nb_parallel,2),"m2")
    print("Depth 2nd floor :",round(best_solution_cma[3],2),"m")
    print("---")
    print("Total volume :",round(Q/best_solution_cma[0]*vdns1.Nb_parallel*best_solution_cma[1]+Q/best_solution_cma[2]*vdns2.Nb_parallel*best_solution_cma[3],2),"m3")
    
    treatment_chain_french_TW_opti = treatment.Treatment_Chain(pathway_french_TW, best_solution_cma, Cin, Cobj, Q)
    output_function_values = treatment_chain_french_TW_opti.Output_Function(best_solution_cma)
    
    print("")
    print("Outlet concentration:")
    print("")
    print("TSS (mgTSS/L):", round(output_function_values[0], 2))
    print("BOD5 (mgO2/L):", round(output_function_values[1], 2))
    print("TKN (mgTKN/L):", round(output_function_values[2], 2))
    print("COD (mgO2/L):", round(output_function_values[3] + output_function_values[4] + output_function_values[5], 2))
        
    print("")
    print("Checking constraints values for best solution found by CMA-ES...")
    print("---")
    output = Cin
    for index, process in enumerate(pathway_french_TW):
        V_values = best_solution_cma[index * 2: (index + 1) * 2]
        
        constraint_TSS_pc = 100 * output[0] * V_values[0] / process.Lim_TSS 
        print("TSS loading contraint stage", index+1, ":", round(constraint_TSS_pc,2), "%")
        constraint_BOD_pc = 100 * output[1] * V_values[0] / process.Lim_BOD
        print("BOD loading contraint stage", index+1, ":", round(constraint_BOD_pc,2), "%")
        constraint_TKN_pc = 100 * output[2] * V_values[0] / process.Lim_TKN
        print("TKN loading contraint stage", index+1, ":", round(constraint_TKN_pc,2), "%")
        constraint_COD_pc = 100 * (output[3]+output[4]+output[5]) * V_values[0] / process.Lim_COD
        print("COD loading contraint stage", index+1, ":", round(constraint_COD_pc,2), "%")
        print("---")
        output = process.Reduction_Function(V_values, output)

    print("Outlet TSS deviation:",round(-(Cobj[0] - output_function_values[0]),2),"mgTSS/L")
    print("Outlet BOD5 deviation:",round(-(Cobj[1] - output_function_values[1]),2),"mgO2/L")
    print("Outlet TKN deviation:",round(-(Cobj[2] - output_function_values[2]),2),"mgTKN/L")
    print("Outlet COD deviation:",round(-(Cobj[3] - (output_function_values[3]+output_function_values[4]+output_function_values[5])),2),"mgO2/L")

################################################################################

if __name__ == '__main__':
    main()
