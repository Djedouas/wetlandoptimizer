import math
import numpy as np
import yaml
import itertools

class Process:
    """
    A class to represent a treatment wetland process. 

    Attributes
    ----------
    Name : str
        Name of the process.
    Lim_TSS : float
        Surfacic TSS load limit (gTSS/m2/day).
    Lim_BOD : float
        Surfacic BOD5 load limit (gO2/m2/day).
    Lim_TKN : float
        Surfacic TKN load limit (gTKN/m2/day).
    Lim_COD : float
        Surfacic COD load limit (gO2/m2/day).
    Nb_parallel : int
        Number of parallel processes for load alternance.
    Material : str
        Material used (gravel or sand).
    Mat_cost : float
        Ton of material cost compared to the cost of one ton of gravel ((Tmat€/m3)/(Tgravel€/m3)).
    Xmin : float
        Hydraulic minimum surface loading rate (m/day).
    Xmax : float
        Hydraulic maximum surface loading rate (m/day).
    Zmin : float
        Minimum unsaturated depth (m).
    Zmax : float
        Maximum unsaturated depth (m).
    Amin :
        Minimum saturated depth (m).
    Amax :
        Maximum saturated depth (m).
    param_TSS : float
        Reduction parameter related to TSS.
    param_BOD : float
        Reduction parameter related to BOD.
    param_TKN_a : float
        Reduction parameter related to TKN (a).
    param_TKN_b : float
        Reduction parameter related to TKN (b).
    param_CODsb : float
        Reduction parameter related to CODsb.
    Cin : list
        Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
    Cobj : list
        Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
    V_values : list
        Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
    Q : float
        Flow rate (m3/day).

    Methods
    -------
    Reduction_Function(V_values, Cin, Q):
        Reduction function for the process.
    Create_Constraint_TSS(V_values, Cin):
        Create TSS constraint, corresponding to the difference between the TSS load constraint and the actual TSS inlet load.
    Create_Constraint_BOD(V_values, Cin):
        Create BOD5 constraint, corresponding to the difference between the DOB5 load constraint and the actual BOD5 inlet load.
    Create_Constraint_TKN(V_values, Cin):
        Create TKN constraint, corresponding to the difference between the TKN load constraint and the actual TKN inlet load.
    Create_Constraint_COD(V_values, Cin):
        Create CODt constraint, corresponding to the difference between the CODt load constraint and the actual CODt inlet load.
    Surface_Area_Function(V_values, Q):
        Calculate the total surface area of the process.  
    Depth_Function_Unsat(V_values):
        Return the unsaturated depth of the process.
    Depth_Function_Sat(self, V_values):
        Return the saturated depth of the process.    
    Depth_Function_Tot(self, V_values):
        Return the total depth of the process.    
    Volume_Function_Unsat(V_values, Q):
        Calculate the unsaturated volume of the process.
    Volume_Function_Sat(self, V_values, Q):
        Calculate the saturated volume of the process.
    Volume_Function_Tot(self, V_values, Q):
        Calculate the total volume of the process.
    Supplementary_Objective_Function(V_values, Cin, Cobj):
        Supplementary objective function for the process.
    Validate_Position(self):
        Defines the rules for the position of the process in the treatment train.
    """
    def __init__(self, Name, Lim_TSS, Lim_BOD, Lim_TKN, Lim_COD, Nb_parallel, Material, Mat_cost, Xmin, Xmax, Zmin, Zmax, Amin, Amax, Param_TSS, Param_BOD, Param_TKN_a, Param_TKN_b, Param_CODsb, Cin=[], Cobj=[],V_values=[], Q=None):
        """
        Initialize a Process object.
        
        Parameters
        ----------
        Name : str
            Name of the process.
        Lim_TSS : float
            Surfacic TSS load limit (gTSS/m2/day).
        Lim_BOD : float
            Surfacic BOD5 load limit (gO2/m2/day).
        Lim_TKN : float
            Surfacic TKN load limit (gTKN/m2/day).
        Lim_COD : float
            Surfacic COD load limit (gO2/m2/day).
        Nb_parallel : int
            Number of parallel processes for load alternance.
        Material : str
            Material used (gravel or sand).
        Mat_cost : float
            Ton of material cost compared to the cost of one ton of gravel ((Tmat€/m3)/(Tgravel€/m3)).
        Xmin : float
            Hydraulic minimum surface loading rate (m/day).
        Xmax : float
            Hydraulic maximum surface loading rate (m/day).
        Zmin : float
            Minimum depth (m).
        Zmax : float
            Maximum depth (m).
        Amin :
            Minimum saturated depth (m).
        Amax :
            Maximum saturated depth (m).
        param_TSS : float
            Reduction parameter related to TSS.
        param_BOD : float
            Reduction parameter related to BOD.
        param_TKN_a : float
            Reduction parameter related to TKN (a).
        param_TKN_b : float
            Reduction parameter related to TKN (b).
        param_CODsb : float
            Reduction parameter related to CODsb.
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Q : float
            Flow rate (m3/day).
        """
        self.Name = Name 
        self.Lim_TSS = Lim_TSS # g/m2/day
        self.Lim_BOD = Lim_BOD # g/m2/day
        self.Lim_TKN = Lim_TKN # g/m2/day
        self.Lim_COD = Lim_COD # g/m2/day
        self.Nb_parallel = Nb_parallel
        self.Material = Material
        self.Mat_cost = Mat_cost # (Tmat€/m3)/(Tgravel€/m3)
        self.Xmin = Xmin # m/day
        self.Xmax = Xmax # m/day
        self.Zmin = Zmin # m
        self.Zmax = Zmax # m
        self.Amax = Amax # %
        self.Amin = Amin  # %
        self.Param_TSS = Param_TSS
        self.Param_BOD = Param_BOD
        self.Param_TKN_a = Param_TKN_a
        self.Param_TKN_b = Param_TKN_b
        self.Param_CODsb = Param_CODsb
        self.Cin = Cin # g/m3
        self.Cobj = Cobj # g/m3
        self.V_values = V_values  
        self.Q = Q # m3/day
    
    def Reduction_Function(self, V_values, Cin, Q):
        """
        Reduction function for the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Q : float
            Flow rate (m3/day).

        Returns
        -------
        list
            Output concentrations after reduction ([TSS]out : Cout[0] (mgTSS/L), [BOD5]out : Cout[1] (mgDOD5/L), [TKN]out : Cout[2] (mgTKN/L), [CODdb]out : Cout[3] (mgCODdb/L), [CODdi]out : Cout[4] (mgCODdi/L), [CODp]out : Cout[5] (mgCODp/L), [NO3]out : Cout[4] (mgNO3/L), [TN]out : Cout[5] (mgN/L)).
        """
        pass
    
    def Create_Constraint_TSS(self, V_values, Cin):
        """
        Create TSS constraint, corresponding to the difference between the TSS load constraint and the actual TSS inlet load.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        
        Returns
        -------
        float
            TSS constraint value.
        """
        return self.Lim_TSS - Cin[0] * V_values[0]
      
    def Create_Constraint_BOD(self, V_values, Cin):
        """
        Create BOD5 constraint, corresponding to the difference between the BOD5 load constraint and the actual BOD5 inlet load.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        
        Returns
        -------
        float
            BOD5 constraint value.
        """
        return self.Lim_BOD - Cin[1] * V_values[0]
    
    def Create_Constraint_TKN(self, V_values, Cin):
        """
        Create TKN constraint, corresponding to the difference between the TKN load constraint and the actual TKN inlet load.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        
        Returns
        -------
        float
            TKN constraint value.
        """
        return self.Lim_TKN - Cin[2] * V_values[0]
    
    def Create_Constraint_COD(self, V_values, Cin):
        """
        Create CODt constraint, corresponding to the difference between the CODt load constraint and the actual CODt inlet load.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        
        Returns
        -------
        float
            CODt constraint value.
        """
        return self.Lim_COD - (Cin[3] + Cin[4] + Cin[5]) * V_values[0]
   
    def Surface_Area_Function(self, V_values, Q):
        """
        Calculate the total surface area of the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Q : float
            Flow rate (m3/day).
        
        Returns
        -------
        float
            Surface area of the process (m2).
        """
        return Q * self.Nb_parallel / V_values[0]
      
    def Depth_Function_Unsat(self, V_values):
        """
        Return the unsaturated depth of the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        
        Returns
        -------
        float
            Unsaturated depth of the process (m).
        """
        return V_values[1]
     
    def Depth_Function_Sat(self, V_values):
        """
        Return the saturated depth of the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        
        Returns
        -------
        float
            Saturated depth of the process (m).
        """
        return V_values[2]
    
    def Depth_Function_Tot(self, V_values):
        """
        Calculate the total depth of the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        
        Returns
        -------
        float
            Total depth of the process (m).
        """
        return self.Depth_Function_Unsat(V_values) + self.Depth_Function_Sat(V_values)
      
    def Volume_Function_Unsat(self, V_values, Q):
        """
        Calculate the unsaturated volume of the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Q : float
            Flow rate (m3/day).
        
        Returns
        -------
        float
            Unsaturated volume of the process (m3).
        """
        return self.Surface_Area_Function(V_values, Q) * self.Depth_Function_Unsat(V_values)

    def Volume_Function_Sat(self, V_values, Q):
        """
        Calculate the saturated volume of the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Q : float
            Flow rate (m3/day).
        
        Returns
        -------
        float
            Saturated volume of the process (m3).
        """
        return self.Surface_Area_Function(V_values, Q) * self.Depth_Function_Sat(V_values)
    
    def Volume_Function_Tot(self, V_values, Q):
        """
        Calculate the total volume of the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Q : float
            Flow rate (m3/day).
        
        Returns
        -------
        float
            Total volume of the process (m3).
        """
        return self.Volume_Function_Unsat(V_values, Q) + self.Volume_Function_Sat(V_values, Q)
    
    def Supplementary_Objective_Function(self, V_values, Cin, Cobj):
        """
        Supplementary objective function for the process.
        
        Parameters
        ----------
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
        
        Returns
        -------
        float
            Value of the supplementary objective function if there is one.
        """
        pass

    def Validate_Position(self):
        """
        Defines the rules for the position of the process in the treatment train.

        Returns
        -------
        bool
            True if the position is valid, otherwise False.
        """
        return True      

################################################################################

class VdNS1(Process):
    """
    A class to represent a first stage type treatment wetland process (VdNS1). 

    Attributes
    ----------
    Name : str
        Name of the process.
    Lim_TSS : float
        Surfacic TSS load limit (gTSS/m2/day).
    Lim_BOD : float
        Surfacic BOD5 load limit (gO2/m2/day).
    Lim_TKN : float
        Surfacic TKN load limit (gTKN/m2/day).
    Lim_COD : float
        Surfacic COD load limit (gO2/m2/day).
    Nb_parallel : int
        Number of parallel processes for load alternance.
    Material : str
        Material used (gravel or sand).
    Mat_cost : float
        Ton of material cost compared to the cost of one ton of gravel ((Tmat€/m3)/(Tgravel€/m3)).
    Xmin : float
        Hydraulic minimum surface loading rate (m/day).
    Xmax : float
        Hydraulic maximum surface loading rate (m/day).
    Zmin : float
        Minimum unsaturated depth (m).
    Zmax : float
        Maximum unsaturated depth (m).
    Amin :
        Minimum saturated depth (m).
    Amax :
        Maximum saturated depth (m).
    param_TSS : float
        Reduction parameter related to TSS.
    param_BOD : float
        Reduction parameter related to BOD.
    param_TKN_a : float
        Reduction parameter related to TKN (a).
    param_TKN_b : float
        Reduction parameter related to TKN (b).
    param_CODsb : float
        Reduction parameter related to CODsb.
    Cin : list
        Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
    Cobj : list
        Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
    V_values : list
        Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
    Q : float
        Flow rate (m3/day).

    Methods
    -------
    Supplementary_Objective_Function(V_values, Cin, Cobj):
        Supplementary objective function for the process.
    Reduction_Function(V_values, Cin, Q):
        Reduction function for the process.
    Create_Constraint_TSS(V_values, Cin):
        Create TSS constraint, corresponding to the difference between the TSS load constraint and the actual TSS inlet load.
    Create_Constraint_BOD(V_values, Cin):
        Create BOD5 constraint, corresponding to the difference between the DOB5 load constraint and the actual BOD5 inlet load.
    Create_Constraint_TKN(V_values, Cin):
        Create TKN constraint, corresponding to the difference between the TKN load constraint and the actual TKN inlet load.
    Create_Constraint_COD(V_values, Cin):
        Create CODt constraint, corresponding to the difference between the CODt load constraint and the actual CODt inlet load.
    Surface_Area_Function(V_values, Q):
        Calculate the total surface area of the process.  
    Depth_Function_Unsat(V_values):
        Return the unsaturated depth of the process.
    Depth_Function_Sat(self, V_values):
        Return the saturated depth of the process.    
    Depth_Function_Tot(self, V_values):
        Return the total depth of the process.    
    Volume_Function_Unsat(V_values, Q):
        Calculate the unsaturated volume of the process.
    Volume_Function_Sat(self, V_values, Q):
        Calculate the saturated volume of the process.
    Volume_Function_Tot(self, V_values, Q):
        Calculate the total volume of the process.
    Optimal_COD_Load(x):
        Calculate the optimal CODt load, as a function of the TKN objective concentration.
    Validate_Position(self):
        Defines the rules for the position of the process in the treatment train.
    """
    def __init__(self, Name, Lim_TSS, Lim_BOD, Lim_TKN, Lim_COD, Nb_parallel, Material, Mat_cost, Xmin, Xmax, Zmin, Zmax, Amin, Amax, Param_TSS, Param_BOD, Param_TKN_a, Param_TKN_b, Param_CODsb, Cin=[], Cobj=[], V_values=[], Q=None):
        """
        Initialize a VdNS1 object.

        This method calls the constructor of the Process class using `super()`.

        Parameters
        ----------
        Name : str
            Name of the process.
        Lim_TSS : float
            Surfacic TSS load limit (gTSS/m2/day).
        Lim_BOD : float
            Surfacic BOD5 load limit (gO2/m2/day).
        Lim_TKN : float
            Surfacic TKN load limit (gTKN/m2/day).
        Lim_COD : float
            Surfacic COD load limit (gO2/m2/day).
        Nb_parallel : int
            Number of parallel processes for load alternance.
        Material : str
            Material used (gravel or sand).
        Mat_cost : float
            Ton of material cost compared to the cost of one ton of gravel ((Tmat€/m3)/(Tgravel€/m3)).
        Xmin : float
            Hydraulic minimum surface loading rate (m/day).
        Xmax : float
            Hydraulic maximum surface loading rate (m/day).
        Zmin : float
            Minimum unsaturated depth (m).
        Zmax : float
            Maximum unsaturated depth (m).
        Amin :
            Minimum saturated depth (m).
        Amax :
            Maximum saturated depth (m).
        param_TSS : float
            Reduction parameter related to TSS.
        param_BOD : float
            Reduction parameter related to BOD.
        param_TKN_a : float
            Reduction parameter related to TKN (a).
        param_TKN_b : float
            Reduction parameter related to TKN (b).
        param_CODsb : float
            Reduction parameter related to CODsb.
        Cin : list
            Input concentrations ([TSS]in : Cin[0] (gTSS/m3), [BOD5]in : Cin[1] (gO2/m3), [TKN]in : Cin[2] (gTKN/m3), [CODdb]in : Cin[3] (gO2/m2), [CODdi]in : Cin[4] (gO2/m3), [CODp]in : Cin[5] (gO2/m3)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (gTSS/m3), [BOD5]obj : Cobj[1] (gO2/m3), [TKN]obj : Cobj[2] (gTKN/m3), [CODt]obj : Cobj[3] (gO2/m2)).
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Q : float
            Flow rate (m3/day).
        """
        super().__init__(Name, Lim_TSS, Lim_BOD, Lim_TKN, Lim_COD, Nb_parallel, Material, Mat_cost, Xmin, Xmax, Zmin, Zmax, Amin, Amax, Param_TSS, Param_BOD, Param_TKN_a, Param_TKN_b, Param_CODsb, Cin, Cobj, V_values, Q)

    def Reduction_Function(self, V_values, Cin, Q) :
        """
        Reduction functions for the first stage type VdNS process.
        
        Parameters
        ----------
        V_values : list
            VdNS1 volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Q : float
            Flow rate (m3/day).

        Returns
        -------
        list
            Output concentrations after reduction ([TSS]out : Cout[0] (mgTSS/L), [BOD5]out : Cout[1] (mgDOD5/L), [TKN]out : Cout[2] (mgTKN/L), [CODdb]out : Cout[3] (mgCODdb/L), [CODdi]out : Cout[4] (mgCODdi/L), [CODp]out : Cout[5] (mgCODp/L), [NO3]out : Cout[4] (mgNO3/L), [TN]out : Cout[5] (mgN/L)).
        """
        TSS_out1 = Cin[0] * self.Param_TSS
        if TSS_out1 < 0:
            TSS_out1 = 0
        
        BOD5_out1 = Cin[1] * self.Param_BOD
        if BOD5_out1 < 0:
            BOD5_out1 = 0
        
        TKN_out1 = Cin[2] - self.Param_TKN_a * ((Cin[2] * V_values[0]) ** self.Param_TKN_b) / V_values[0]
        if TKN_out1 < 0:
            TKN_out1 = 0
        if TKN_out1 > Cin[2]:
            TKN_out1 = Cin[2]
    
        CODsb_out1 = Cin[3] * math.exp(-self.Param_CODsb * V_values[1])
        if CODsb_out1 < 0:
            CODsb_out1 = 0
    
        CODsi_out1 = Cin[4]
    
        CODp_out1 = Cin[5] * self.Param_TSS
        if CODp_out1 < 0:
            CODp_out1 = 0
    
        if CODsb_out1 + CODsi_out1 + CODp_out1 > Cin[3] + Cin[4] + Cin[5]:
            CODsb_out1 = Cin[3]
            CODsi_out1 = Cin[4]
            CODp_out1 = Cin[5]

        NO3_out1 = Cin[6] + (Cin[2] - TKN_out1)
        if NO3_out1 < 0:
            NO3_out1 = 0
        
        Cout1 = [TSS_out1, BOD5_out1, TKN_out1, CODsb_out1, CODsi_out1, CODp_out1, NO3_out1]
    
        return Cout1
      
    def Optimal_COD_Load(self, x):
        """
        Calculate the optimal CODt load, as a function of the TKN objective concentration.
        
        Parameters
        ----------
        x : float
            TKN objective concentration (gTKN/m3).
        
        Returns
        -------
        float
            Optimal CODt load.
        """
        if x == None :
            return 350
        else :
            y = (175/6) * x
            if y >= 350:
                return 350
            else:
                return y
          
    def Supplementary_Objective_Function(self, V_values, Cin, Cobj):
        """
        Supplementary objective function for the first stage type VdNS process, corresponding to the normalized squared error between the optimal CODt load and the actual CODt load.
        
        Parameters
        ----------
        V_values : list
            VdNS1 volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
        
        Returns
        -------
        float
            Result of the supplementary objective function.
        """
        optimal_load_COD_1 = self.Optimal_COD_Load(Cobj[2])
        load_COD_1 = (Cin[3] + Cin[4] + Cin[5]) * V_values[0]
        clogg_COD_1 = ((optimal_load_COD_1 - load_COD_1) ** 2) / (optimal_load_COD_1 ** 2)
        return clogg_COD_1 * 200
    
    def Validate_Position(self,combination):
        """
        Defines the rules for the position of the first stage type VdNS process in the treatment train.

        Returns
        -------
        bool
            True if the position is valid, otherwise False.
        """
        for i, process in enumerate(combination):
            if process == self:
                if i != 0:
                    return False
        return True
             
################################################################################

class VdNS2(Process):
    """
    A class to represent a second stage type treatment wetland process (VdNS2). 

    Attributes
    ----------
    Name : str
        Name of the process.
    Lim_TSS : float
        Surfacic TSS load limit (gTSS/m2/day).
    Lim_BOD : float
        Surfacic BOD5 load limit (gO2/m2/day).
    Lim_TKN : float
        Surfacic TKN load limit (gTKN/m2/day).
    Lim_COD : float
        Surfacic COD load limit (gO2/m2/day).
    Nb_parallel : int
        Number of parallel processes for load alternance.
    Material : str
        Material used (gravel or sand).
    Mat_cost : float
        Ton of material cost compared to the cost of one ton of gravel ((Tmat€/m3)/(Tgravel€/m3)).
    Xmin : float
        Hydraulic minimum surface loading rate (m/day).
    Xmax : float
        Hydraulic maximum surface loading rate (m/day).
    Zmin : float
        Minimum unsaturated depth (m).
    Zmax : float
        Maximum unsaturated depth (m).
    Amin :
        Minimum saturated depth (m).
    Amax :
        Maximum saturated depth (m).
    param_TSS : float
        Reduction parameter related to TSS.
    param_BOD : float
        Reduction parameter related to BOD.
    param_TKN_a : float
        Reduction parameter related to TKN (a).
    param_TKN_b : float
        Reduction parameter related to TKN (b).
    param_CODsb : float
        Reduction parameter related to CODsb.
    Cin : list
        Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
    Cobj : list
        Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
    V_values : list
        Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
    Q : float
        Flow rate (m3/day).

    Methods
    -------
    Supplementary_Objective_Function(V_values, Cin, Cobj):
        Supplementary objective function for the process.
    Reduction_Function(V_values, Cin, Q):
        Reduction function for the process.
    Create_Constraint_TSS(V_values, Cin):
        Create TSS constraint, corresponding to the difference between the TSS load constraint and the actual TSS inlet load.
    Create_Constraint_BOD(V_values, Cin):
        Create BOD5 constraint, corresponding to the difference between the DOB5 load constraint and the actual BOD5 inlet load.
    Create_Constraint_TKN(V_values, Cin):
        Create TKN constraint, corresponding to the difference between the TKN load constraint and the actual TKN inlet load.
    Create_Constraint_COD(V_values, Cin):
        Create CODt constraint, corresponding to the difference between the CODt load constraint and the actual CODt inlet load.
    Surface_Area_Function(V_values, Q):
        Calculate the total surface area of the process.  
    Depth_Function_Unsat(V_values):
        Return the unsaturated depth of the process.
    Depth_Function_Sat(self, V_values):
        Return the saturated depth of the process.    
    Depth_Function_Tot(self, V_values):
        Return the total depth of the process.    
    Volume_Function_Unsat(V_values, Q):
        Calculate the unsaturated volume of the process.
    Volume_Function_Sat(self, V_values, Q):
        Calculate the saturated volume of the process.
    Volume_Function_Tot(self, V_values, Q):
        Calculate the total volume of the process.
    Optimal_COD_Load(x):
        Calculate the optimal CODt load, as a function of the TKN objective concentration.
    Validate_Position(self):
        Defines the rules for the position of the process in the treatment train.
    """
    def __init__(self, Name, Lim_TSS, Lim_BOD, Lim_TKN, Lim_COD, Nb_parallel, Material, Mat_cost, Xmin, Xmax, Zmin, Zmax, Amin, Amax, Param_TSS, Param_BOD, Param_TKN_a, Param_TKN_b, Param_CODsb,  Cin=[], Cobj=[], V_values=[], Q=None):
        """
        Initialize a VdNS2 object.

        This method calls the constructor of the Process class using `super()`.

        Parameters
        ----------
        Name : str
            Name of the process.
        Lim_TSS : float
            Surfacic TSS load limit (gTSS/m2/day).
        Lim_BOD : float
            Surfacic BOD5 load limit (gO2/m2/day).
        Lim_TKN : float
            Surfacic TKN load limit (gTKN/m2/day).
        Lim_COD : float
            Surfacic COD load limit (gO2/m2/day).
        Nb_parallel : int
            Number of parallel processes for load alternance.
        Material : str
            Material used (gravel or sand).
        Mat_cost : float
            Ton of material cost compared to the cost of one ton of gravel ((Tmat€/m3)/(Tgravel€/m3)).
        Xmin : float
            Hydraulic minimum surface loading rate (m/day).
        Xmax : float
            Hydraulic maximum surface loading rate (m/day).
        Zmin : float
            Minimum unsaturated depth (m).
        Zmax : float
            Maximum unsaturated depth (m).
        Amin :
            Minimum saturated depth (m).
        Amax :
            Maximum saturated depth (m).
        param_TSS : float
            Reduction parameter related to TSS.
        param_BOD : float
            Reduction parameter related to BOD.
        param_TKN_a : float
            Reduction parameter related to TKN (a).
        param_TKN_b : float
            Reduction parameter related to TKN (b).
        param_CODsb : float
            Reduction parameter related to CODsb.
        Cin : list
            Input concentrations ([TSS]in : Cin[0] (gTSS/m3), [BOD5]in : Cin[1] (gO2/m3), [TKN]in : Cin[2] (gTKN/m3), [CODdb]in : Cin[3] (gO2/m2), [CODdi]in : Cin[4] (gO2/m3), [CODp]in : Cin[5] (gO2/m3)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (gTSS/m3), [BOD5]obj : Cobj[1] (gO2/m3), [TKN]obj : Cobj[2] (gTKN/m3), [CODt]obj : Cobj[3] (gO2/m2)).
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Q : float
            Flow rate (m3/day).
        """
        super().__init__(Name, Lim_TSS, Lim_BOD, Lim_TKN, Lim_COD, Nb_parallel, Material, Mat_cost, Xmin, Xmax, Zmin, Zmax, Amin, Amax, Param_TSS, Param_BOD, Param_TKN_a, Param_TKN_b, Param_CODsb, Cin, Cobj, V_values, Q)

    def Reduction_Function(self, V_values, Cin, Q) :
        """
        Reduction functions for the second stage type VdNS process.
        
        Parameters
        ----------
        V_values : list
            VdNS1 volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Q : float
            Flow rate (m3/day).

        Returns
        -------
        list
            Output concentrations after reduction ([TSS]out : Cout[0] (mgTSS/L), [BOD5]out : Cout[1] (mgDOD5/L), [TKN]out : Cout[2] (mgTKN/L), [CODdb]out : Cout[3] (mgCODdb/L), [CODdi]out : Cout[4] (mgCODdi/L), [CODp]out : Cout[5] (mgCODp/L), [NO3]out : Cout[4] (mgNO3/L), [TN]out : Cout[5] (mgN/L)).
        """
        TSS_out2 = Cin[0] * self.Param_TSS
        if TSS_out2 < 0 :
            TSS_out2 = 0
        
        BOD5_out2 = Cin[1] * self.Param_BOD
        if BOD5_out2 < 0 :
            BOD5_out2 = 0
      
        TKN_out2 = Cin[2] - self.Param_TKN_a * ((Cin[2] * V_values[0]) ** self.Param_TKN_b) / V_values[0]
        if TKN_out2 < 0 :
            TKN_out2 = 0
        if TKN_out2 > Cin[2] :
            TKN_out2 = Cin[2]  
      
        CODsb_out2 = Cin[3] * math.exp(-self.Param_CODsb * V_values[1])
        if CODsb_out2 < 0 :
            CODsb_out2 = 0
      
        CODsi_out2 = Cin[4]
        if CODsi_out2 < 0 :
            CODsi_out2 = 0
        
        CODp_out2 = Cin[5] * self.Param_TSS
        if CODp_out2 < 0 :
            CODp_out2 = 0
      
        if CODsb_out2 + CODsi_out2 + CODp_out2 > Cin[3] +  Cin[4] +  Cin[5] :
            CODsb_out2 = Cin[3]
            CODsi_out2 = Cin[4]
            CODp_out2 = Cin[5]
        
        NO3_out2 = Cin[6] + (Cin[2] - TKN_out2)
        if NO3_out2 < 0:
            NO3_out2 = 0

        Cout2 = [TSS_out2,BOD5_out2,TKN_out2,CODsb_out2,CODsi_out2,CODp_out2, NO3_out2]
        
        return Cout2
      
    def Supplementary_Objective_Function(self, V_values, Cin, Cobj):
        """
        Supplementary objective function for the second stage type VdNS process, which is equal to zero for this process.
        
        Parameters
        ----------
        V_values : list
            VdNS1 volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
        
        Returns
        -------
        float
            Result of the supplementary objective function.
        """
        return 0
    
    def Validate_Position(self, combination):
        """
        Defines the rules for the position of the second stage type VdNS process in the treatment train.

        Returns
        -------
        bool
            True if the position is valid, otherwise False.
        """
        for i, process in enumerate(combination):
            if process == self:
                if i == 0:
                    return False
                return True
        return True
    
################################################################################

class VdNSS(Process):
    """
    A class to represent a first stage with a saturated bottom type treatment wetland process (VdNSS). 

    Attributes
    ----------
    Name : str
        Name of the process.
    Lim_TSS : float
        Surfacic TSS load limit (gTSS/m2/day).
    Lim_BOD : float
        Surfacic BOD5 load limit (gO2/m2/day).
    Lim_TKN : float
        Surfacic TKN load limit (gTKN/m2/day).
    Lim_COD : float
        Surfacic COD load limit (gO2/m2/day).
    Nb_parallel : int
        Number of parallel processes for load alternance.
    Material : str
        Material used (gravel or sand).
    Mat_cost : float
        Ton of material cost compared to the cost of one ton of gravel ((Tmat€/m3)/(Tgravel€/m3)).
    Xmin : float
        Hydraulic minimum surface loading rate (m/day).
    Xmax : float
        Hydraulic maximum surface loading rate (m/day).
    Zmin : float
        Minimum unsaturated depth (m).
    Zmax : float
        Maximum unsaturated depth (m).
    Amin :
        Minimum saturated depth (m).
    Amax :
        Maximum saturated depth (m).
    param_TSS : float
        Reduction parameter related to TSS.
    param_BOD : float
        Reduction parameter related to BOD.
    param_TKN_a : float
        Reduction parameter related to TKN (a).
    param_TKN_b : float
        Reduction parameter related to TKN (b).
    param_CODsb : float
        Reduction parameter related to CODsb.
    Cin : list
        Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
    Cobj : list
        Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
    V_values : list
        Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
    Q : float
        Flow rate (m3/day).

    Methods
    -------
    Supplementary_Objective_Function(V_values, Cin, Cobj):
        Supplementary objective function for the process.
    Reduction_Function(V_values, Cin, Q):
        Reduction function for the process.
    Create_Constraint_TSS(V_values, Cin):
        Create TSS constraint, corresponding to the difference between the TSS load constraint and the actual TSS inlet load.
    Create_Constraint_BOD(V_values, Cin):
        Create BOD5 constraint, corresponding to the difference between the DOB5 load constraint and the actual BOD5 inlet load.
    Create_Constraint_TKN(V_values, Cin):
        Create TKN constraint, corresponding to the difference between the TKN load constraint and the actual TKN inlet load.
    Create_Constraint_COD(V_values, Cin):
        Create CODt constraint, corresponding to the difference between the CODt load constraint and the actual CODt inlet load.
    Surface_Area_Function(V_values, Q):
        Calculate the total surface area of the process.  
    Depth_Function_Unsat(V_values):
        Return the unsaturated depth of the process.
    Depth_Function_Sat(self, V_values):
        Return the saturated depth of the process.    
    Depth_Function_Tot(self, V_values):
        Return the total depth of the process.    
    Volume_Function_Unsat(V_values, Q):
        Calculate the unsaturated volume of the process.
    Volume_Function_Sat(self, V_values, Q):
        Calculate the saturated volume of the process.
    Volume_Function_Tot(self, V_values, Q):
        Calculate the total volume of the process.
    Validate_Position(self):
        Defines the rules for the position of the process in the treatment train.
    """ 
    def __init__(self, Name, Lim_TSS, Lim_BOD, Lim_TKN, Lim_COD, Nb_parallel, Material, Mat_cost, Xmin, Xmax, Zmin, Zmax, Amin, Amax, Param_TSS, Param_BOD, Param_TKN_a, Param_TKN_b, Param_CODsb, Cin=[], Cobj=[], V_values=[], Q=None):
        """
        Initialize a VdNSS object.

        This method calls the constructor of the Process class using `super()`.

        Parameters
        ----------
        Name : str
            Name of the process.
        Lim_TSS : float
            Surfacic TSS load limit (gTSS/m2/day).
        Lim_BOD : float
            Surfacic BOD5 load limit (gO2/m2/day).
        Lim_TKN : float
            Surfacic TKN load limit (gTKN/m2/day).
        Lim_COD : float
            Surfacic COD load limit (gO2/m2/day).
        Nb_parallel : int
            Number of parallel processes for load alternance.
        Material : str
            Material used (gravel or sand).
        Mat_cost : float
            Ton of material cost compared to the cost of one ton of gravel ((Tmat€/m3)/(Tgravel€/m3)).
        Xmin : float
            Hydraulic minimum surface loading rate (m/day).
        Xmax : float
            Hydraulic maximum surface loading rate (m/day).
        Zmin : float
            Minimum unsaturated depth (m).
        Zmax : float
            Maximum unsaturated depth (m).
        Amin :
            Minimum saturated depth (m).
        Amax :
            Maximum saturated depth (m).
        param_TSS : float
            Reduction parameter related to TSS.
        param_BOD : float
            Reduction parameter related to BOD.
        param_TKN_a : float
            Reduction parameter related to TKN (a).
        param_TKN_b : float
            Reduction parameter related to TKN (b).
        param_CODsb : float
            Reduction parameter related to CODsb.
        Cin : list
            Input concentrations ([TSS]in : Cin[0] (gTSS/m3), [BOD5]in : Cin[1] (gO2/m3), [TKN]in : Cin[2] (gTKN/m3), [CODdb]in : Cin[3] (gO2/m2), [CODdi]in : Cin[4] (gO2/m3), [CODp]in : Cin[5] (gO2/m3)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (gTSS/m3), [BOD5]obj : Cobj[1] (gO2/m3), [TKN]obj : Cobj[2] (gTKN/m3), [CODt]obj : Cobj[3] (gO2/m2)).
        V_values : list
            Process volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Q : float
            Flow rate (m3/day).
        """
        super().__init__(Name, Lim_TSS, Lim_BOD, Lim_TKN, Lim_COD, Nb_parallel, Material, Mat_cost, Xmin, Xmax, Zmin, Zmax, Amin, Amax, Param_TSS, Param_BOD, Param_TKN_a, Param_TKN_b, Param_CODsb, Cin, Cobj, V_values, Q)

    def Reduction_Function(self, V_values, Cin, Q) :
        """
        Reduction functions for the first stage with a saturated bottom type process.
        
        Parameters
        ----------
        V_values : list
            VdNS1 volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Q : float
            Flow rate (m3/day).

        Returns
        -------
        list
            Output concentrations after reduction ([TSS]out : Cout[0] (mgTSS/L), [BOD5]out : Cout[1] (mgDOD5/L), [TKN]out : Cout[2] (mgTKN/L), [CODdb]out : Cout[3] (mgCODdb/L), [CODdi]out : Cout[4] (mgCODdi/L), [CODp]out : Cout[5] (mgCODp/L), [NO3]out : Cout[4] (mgNO3/L), [TN]out : Cout[5] (mgN/L)).
        """
        TSS_out1 = Cin[0] * self.Param_TSS
        if TSS_out1 < 0:
            TSS_out1 = 0
        
        BOD5_out1 = Cin[1] * self.Param_BOD
        if BOD5_out1 < 0:
            BOD5_out1 = 0
        
        TKN_out1 = Cin[2] - self.Param_TKN_a * ((Cin[2] * V_values[0]) ** self.Param_TKN_b) / V_values[0]
        if TKN_out1 < 0:
            TKN_out1 = 0
        if TKN_out1 > Cin[2]:
            TKN_out1 = Cin[2]
    
        CODsb_out_unsat1 = Cin[3] * math.exp(-self.Param_CODsb * V_values[1])
        CODsb_out1 = CODsb_out_unsat1 * math.exp(-self.Param_CODsb * V_values[2]) 
        if CODsb_out1 < 0:
            CODsb_out1 = 0

        CODsi_out1 = Cin[4]
    
        CODp_out1 = Cin[5] * self.Param_TSS
        if CODp_out1 < 0:
            CODp_out1 = 0
    
        if CODsb_out1 + CODsi_out1 + CODp_out1 > Cin[3] + Cin[4] + Cin[5]:
            CODsb_out1 = Cin[3]
            CODsi_out1 = Cin[4]
            CODp_out1 = Cin[5]

        ratio = (CODsb_out_unsat1 + CODsi_out1 + CODp_out1) / (Cin[6] + (Cin[2] - TKN_out1))
 
        NO3_out1 = (Cin[6] + (Cin[2] - TKN_out1)) / ((1 + ((1.98*math.exp(0.49 * ratio)-1) * V_values[2] * Q/V_values[0] * 0.3) / ( Q * (1.33 *  V_values[2]+0.68)))**(1.33 *  V_values[2] + 0.68))
        if NO3_out1 < 0:
            NO3_out1 = 0

        Cout1 = [TSS_out1, BOD5_out1, TKN_out1, CODsb_out1, CODsi_out1, CODp_out1, NO3_out1]
    
        return Cout1
          
    def Supplementary_Objective_Function(self, V_values, Cin, Cobj):
        """
        Supplementary objective function for the first stage with a saturated bottom type process, corresponding to the normalized squared error between the optimal CODt load and the actual CODt load.
        
        Parameters
        ----------
        V_values : list
            VdNS1 volume values (Q / surface area : V[0] (m/day) and depth : V[1] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
        
        Returns
        -------
        float
            Result of the supplementary objective function.
        """ 
        return 0

    def Validate_Position(self,combination) :
        """
        Defines the rules for the position of the first stage with a staurated bottom type process in the treatment train.

        Returns
        -------
        bool
            True if the position is valid, otherwise False.
        """        
        for i, process in enumerate(combination):
            if process == self:
                if i != 0:
                    return False
        return True

################################################################################

class Treatment_Train:
    """
    A class to represent a two stages wastewater treatment train.

    Attributes
    ----------
    pathway : list
        List of treatment processes in the train.
    Cin : list
        Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
    Cobj : list
        Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
    Q : float
        Flow rate (m3/day).
    V : list
        Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).

    Methods
    -------
    Output_Function(V, Q):
        Calculate the total treatment train output concentration of pollutants, by using each of the two stages reduction functions.
    Create_Constraints_TSS(V, Q):
        Combinate total treatment train load constraints for TSS.
    Create_Constraints_BOD(V, Q):
        Combinate total treatment train load constraints for BOD5.
    Create_Constraints_TKN(V, Q):
        Combinate total treatment train load constraints for TKN.
    Create_Constraints_COD(V, Q):
        Combinate total treatment train load constraints for CODt.
    Create_Constraints_TSSout(V, Cobj, Q):
        Create total treatment train output constraints for TSS, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.
    Create_Constraints_BODout(V, Cobj, Q):
        Create total treatment train output constraints for BOD5, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.
    Create_Constraints_TKNout(V, Cobj, Q):
        Create total treatment train output constraints for TKN, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.
    Create_Constraints_CODout(V, Cobj, Q):
        Create total treatment train output constraints for CODt, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.
    Create_Constraints_NO3out(V, Cobj, Q):
        Create total treatment train output constraints for NO3, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.
    Create_Constraints_TNout(V, Cobj, Q):
        Create total treatment train output constraints for TN, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.
    Total_Volume_Function(V, Q):
        Calculate the total volume of the treatment train.
    Total_Surface_Area_Function(self, V, Q) :
        Calculate the total surrface area of the treatment train.
    Create_Objective_Function(V, Cin, Cobj, Q):
        Combinate the different objective function subparts (volume with economic weighting and supplementary objective function for the first stage) for optimization.
    Single_Objective_Function(V, Cin, Cobj, Q):
        Sum of the different results of the objective function subparts and constraints, to optimize a single fitness value.
    """
    def __init__(self, pathway, V, Cin, Cobj, Q):
        """
        Initialize the treatment train with given parameters.

        Parameters
        ----------
        pathway : list
            List of treatment processes in the treatment train (VdNS1: pathway[0], VdNS2: pathway[1]).
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
        Q : float
            Flow rate (m3/day).
        """
        self.pathway = pathway
        self.Cin = Cin # g/m3
        self.Cobj = Cobj # g/m3
        self.Q = Q # m3/day
        self.V = V # m3

    def Output_Function(self, V, Q):
        """
        Calculate the total treatment train output concentration of pollutants, by using each of the two stages reduction functions.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).

        Returns
        -------
        output : list
            Final output concentration of pollutants ([TSS]out : Cout[0] (mgTSS/L), [BOD5]out : Cout[1] (mgDOD5/L), [TKN]out : Cout[2] (mgTKN/L), [CODdb]out : Cout[3] (mgCODdb/L), [CODdi]out : Cout[4] (mgCODdi/L), [CODp]out : Cout[5] (mgCODp/L), [NO3]out : Cout[4] (mgNO3/L), [TN]out : Cout[5] (mgN/L)).
        """
        output = self.Cin
        for index, process in enumerate(self.pathway):
            output = process.Reduction_Function(V[index * 3: (index + 1) * 3], output, Q)
        return output
      
    def Create_Constraints_TSS(self, V, Q) :
        """
        Combinate total treatment train load constraints for TSS.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).

        Returns
        -------
        constraints : list
            List of TSS load constraint values.
        """
        constraints = []
        output = self.Cin
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            constraint = process.Create_Constraint_TSS(V_values, output)
            constraints.append(constraint)
            output = process.Reduction_Function(V_values, output, Q)
        return constraints
      
    def Create_Constraints_BOD(self, V, Q) :
        """
        Combinate total treatment train load constraints for BOD5.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).

        Returns
        -------
        constraints : list
            List of BOD5 load constraint values.
        """
        constraints = []
        output = self.Cin
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            constraint = process.Create_Constraint_BOD(V_values, output)
            constraints.append(constraint)
            output = process.Reduction_Function(V_values, output, Q)
        return constraints

    def Create_Constraints_TKN(self, V, Q) :
        """
        Combinate total treatment train load constraints for TKN.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).

        Returns
        -------
        constraints : list
            List of TKN load constraint values.
        """
        constraints = []
        output = self.Cin
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            constraint = process.Create_Constraint_TKN(V_values, output)
            constraints.append(constraint)
            output = process.Reduction_Function(V_values, output, Q)
        return constraints
      
    def Create_Constraints_COD(self, V, Q) :
        """
        Combinate total treatment train load constraints for CODt.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).

        Returns
        -------
        constraints : list
            List of CODt load constraint values.
        """
        constraints = []
        output = self.Cin
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            constraint = process.Create_Constraint_COD(V_values, output)
            constraints.append(constraint)
            output = process.Reduction_Function(V_values, output, Q)
        return constraints
    
    def Create_Constraints_TSSout(self, V, Cobj, Q) :
        """
        Create total treatment train output constraints for TSS, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).

        Returns
        -------
        constraint : list
            TSS output constraint value.
        """
        output = self.Cin
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            output = process.Reduction_Function(V_values, output, Q)
        constraint = Cobj[0] - output[0]
        return [constraint]
    
    def Create_Constraints_BODout(self, V, Cobj, Q) :
        """
        Create total treatment train output constraints for BOD5, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).

        Returns
        -------
        constraint : list
            BOD5 output constraint value.
        """
        output = self.Cin
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            output = process.Reduction_Function(V_values, output, Q)
        constraint = Cobj[1] - output[1]
        return [constraint]

    def Create_Constraints_TKNout(self, V, Cobj, Q) :
        """
        Create total treatment train output constraints for TKN, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).

        Returns
        -------
        constraint : list
            TKN output constraint value.
        """
        if Cobj[2] == None :
            return [0]
        else :
            output = self.Cin
            for index, process in enumerate(self.pathway):
                V_values = V[index * 3: (index + 1) * 3]
                output = process.Reduction_Function(V_values, output, Q)
            constraint = Cobj[2] - output[2]
            return [constraint]
      
    def Create_Constraints_CODout(self, V, Cobj, Q) :
        """
        Create total treatment train output constraints for COD, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).

        Returns
        -------
        constraint : list
            COD output constraint value.
        """
        output = self.Cin
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            output = process.Reduction_Function(V_values, output, Q)
        constraint = Cobj[3] - (output[3] + output[4] + output[5])
        return [constraint]
    
    def Create_Constraints_NO3out(self, V, Cobj, Q) :
        """
        Create total treatment train output constraints for NO3, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).

        Returns
        -------
        constraint : list
            NO3 output constraint value.
        """
        if Cobj[4] == None :
            return [0]
        else :
            output = self.Cin
            for index, process in enumerate(self.pathway):
                V_values = V[index * 3: (index + 1) * 3]
                output = process.Reduction_Function(V_values, output, Q)
            constraint = Cobj[4] - output[6]
            return [constraint]

    def Create_Constraints_TNout(self, V, Cobj, Q) :
        """
        Create total treatment train output constraints for TN, corresponding to the difference between the objective outlet concentration and the actuel outlet concentration.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).

        Returns
        -------
        constraint : list
            TN output constraint value.
        """
        if Cobj[5] == None :
            return [0]
        else :    
            output = self.Cin
            for index, process in enumerate(self.pathway):
                V_values = V[index * 3: (index + 1) * 3]
                output = process.Reduction_Function(V_values, output, Q)
            constraint = Cobj[5] - (output[2]+output[6])
            return [constraint]    
    
    def Total_Volume_Function(self, V, Q) :
        """
        Calculate the total volume of the treatment train.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Q : float
            Flow rate (m3/day).

        Returns
        -------
        volume : float
            Total volume of the treatment train (m3).
        """
        volume = 0
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            vol = process.Volume_Function_Tot(V_values,Q)
            volume += vol
        return volume
    
    def Total_Surface_Area_Function(self, V, Q) :
        """
        Calculate the total surface area of the treatment train.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Q : float
            Flow rate (m3/day).

        Returns
        -------
        surface_area : float
            Total surface area of the treatment train (m2).
        """
        surface_area = 0
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            s_a = process.Surface_Area_Function(V_values, Q)
            surface_area += s_a
        return surface_area
     
    def Create_Objective_Function(self, V, Cin, Cobj, Q):
        """
        Combinate the different objective function subparts (volume with economic weighting and supplementary objective function for the first stage) for optimization.

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
        Q : float
            Flow rate (m3/day).

        Returns
        -------
        objective : list
            List of objective function subparts values for optimization (volume with economic weighting for VdNS1: objective[0], supplementary objective function for VdNS1: objective[1], volume with economic weighting for VdNS2[2], supplementary objective function for VdNS2 (null): objective[3]).
        """
        objective = []
        output = self.Cin
        for index, process in enumerate(self.pathway):
            V_values = V[index * 3: (index + 1) * 3]
            vol = process.Volume_Function_Unsat(V_values,Q) * process.Mat_cost + process.Volume_Function_Sat(V_values,Q) * 0.66
            #demander pq 0.66 !!
            objective.append(vol)
            supp = process.Supplementary_Objective_Function(V_values, output, Cobj)
            objective.append(supp)
            output = process.Reduction_Function(V_values, output, Q)
        return objective
      
    def Single_Objective_Function(self, V, Cin, Cobj, Q):
        """
        Sum of the different results of the objective function subparts and constraints, to optimize a single fitness value.  

        Parameters
        ----------
        V : list
            Total treatment train volume values (Q / surface area of the first stage: V[0] (m/day) and depth of the first stage: V[1] (m), Q / surface area of the second stage: V[2] (m/day) and depth of the second stage: V[3] (m)).
        Cin : list
            Input concentrations ([TSS]in1 : Cin[0] (mgTSS/L), [BOD5]in1 : Cin[1] (mgBOD5/L), [TKN]in1 : Cin[2] (mgTKN/L), [CODdb]in1 : Cin[3] (mgCODdb/L), [CODdi]in1 : Cin[4] (mgCODdi/m3), [CODp]in1 : Cin[5] (mgCODp/L), [NO3]in1 : Cin[6] (mgNO3/L)).
        Cobj : list
            Objective concentrations ([TSS]obj : Cobj[0] (mgTSS/L), [BOD5]obj : Cobj[1] (mgBOD5/L), [TKN]obj : Cobj[2] (mgTKN/L), [CODt]obj : Cobj[3] (mgCODt/m3), [NO3]obj : Cobj[4] (mgNO3/L), [TN]obj : Cobj[5] (mgN/L)).
        Q : float
            Flow rate (m3/day).

        Returns
        -------
        fitness : float
            Fitness value.
        """
        constraints_TSS = self.Create_Constraints_TSS(V, Q)
        constraints_BOD = self.Create_Constraints_BOD(V, Q)
        constraints_TKN = self.Create_Constraints_TKN(V, Q)
        constraints_COD = self.Create_Constraints_COD(V, Q)
        constraint_TSSout = self.Create_Constraints_TSSout(V, Cobj, Q)
        constraint_BODout = self.Create_Constraints_BODout(V, Cobj, Q)
        constraint_TKNout = self.Create_Constraints_TKNout(V, Cobj, Q)
        constraint_CODout = self.Create_Constraints_CODout(V, Cobj, Q)
        constraint_NO3out = self.Create_Constraints_NO3out(V, Cobj, Q)
        constraint_TNout = self.Create_Constraints_TNout(V, Cobj, Q)
        constraints = []
        constraints.extend(constraints_TSS)
        constraints.extend(constraints_BOD)
        constraints.extend(constraints_TKN)
        constraints.extend(constraints_COD)
        constraints.extend(constraint_TSSout)
        constraints.extend(constraint_BODout)
        constraints.extend(constraint_TKNout)
        constraints.extend(constraint_CODout)
        constraints.extend(constraint_NO3out)
        constraints.extend(constraint_TNout)
        objectives = self.Create_Objective_Function(V, Cin, Cobj, Q)
        constraint_penalty = 0
        obj_value = 0
        fitness = 0
        for constraint in constraints : 
          if constraint < 0 :
            constraint_penalty += abs(constraint)*10000
        for objectif in objectives :
          obj_value += objectif 
        fitness = 10000 * constraint_penalty + obj_value
        return fitness
    
##################################################################################

class Pathway:
    """
    A class handling the configuration and generation of valid treatment train pathways.

    Attributes
    ----------
    stages_max : float
        Maximum value for the number of stages in series.
    files_max : float
        Maximum value for the number of files in parallel.

    Methods
    -------
    Get_Subclasses(cls):
        Recursively retrieves all subclasses of a given class.
    Possible_Combinations():
        Generates all valid process combinations with exactly stages_max steps.
    """
    with open("config.yaml", "r") as file:
        config = yaml.safe_load(file)

    def __init__(self, stages_max, files_max):
        """
        Initialize the pathway class with given parameters.

        Parameters
        ----------
        stages_max : float
            Maximum value for the number of stages in series.
        files_max : float
            Maximum value for the number of files in parallel.
        """
        self.stages_max = stages_max
        self.files_max = files_max

        subclasses = self.Get_Subclasses(Process)

        self.processes = [
            subclass(**self.config[subclass.__name__])
            for subclass in subclasses
            if subclass.__name__ in self.config
        ]

    @staticmethod

    def Get_Subclasses(cls):
        """
        Recursively retrieves all subclasses of a given class.

        Parameters
        ----------
        cls : type
            The class whose subclasses are to be retrieved.

        Returns
        -------
        list
            List of all found subclasses.
        """
        subclasses = []
        for subclass in cls.__subclasses__():
            subclasses.append(subclass)
            subclasses.extend(Pathway.Get_Subclasses(subclass))
        return subclasses

    def Possible_Combinations(self):
        """
        Generates all valid process combinations with exactly stages_max steps.

        Returns
        -------
        set
            Set of valid process combinations as tuples.
        """
        combinations = set()
        raw_combinations = itertools.product(self.processes, repeat=self.stages_max)
        for combination in raw_combinations:
            # Valider la combinaison en appelant Validate_Position pour chaque procédé
            if all(process.Validate_Position(combination) for process in combination):
                combinations.add(combination)
        return combinations

