import numpy as np

# Import here `odeint` and `interp1d' to improve performance (multi-call)
from mdl import mixed_culture_ode
from scipy.integrate import odeint
from scipy.interpolate import interp1d

def _calculate_carbon(C_range, A_range, dr, MIC):
    """
        Calculate relative carbon content at the MIC. The method calculates
        the `A_range' indexes immediately above and below the MIC, and applies
        them to `C_range'.
        
        Arguments:
            A_range: array. Range of drug concentrations.
            C_range: array. Total carbon uptaken across `A_range'.
            dr: array. Microbial growth across `A_range' (dose-response).
            MIC: float. Minimum Inhibitory Concentration.
    """
    Rel_CiC = C_range/dr
    if np.size(np.where(A_range <= MIC)) == 0:
        L_ID = np.nan
    else:
        L_ID = np.where(A_range <= MIC)[0][-1]
    
    if np.size(np.where(A_range >= MIC)) == 0:
        U_ID = np.nan
    else:
        U_ID = np.where(A_range >= MIC)[0][0]
    
    if np.isnan(L_ID) or np.isnan(U_ID):
        # No MIC defined, return NaN
        return np.nan
    else:
        L_CiC = Rel_CiC[L_ID]
        U_Cic = Rel_CiC[U_ID]
    return np.mean([L_CiC, U_Cic])

def ode_solver(u0, t, p, Competitor_Type='Resistant'):
    """
        Solves ODE system in `f' and returns the solution.
        
        Arguments:
            f: Function containing an ODE system.
            u0: List containing the initial guesses for all the _variables_
                  of the model.
            t: Timespan to solve the ODE system as an array.
            p: List containing all the _parameters_ of the model.
    """
    solution, _ = odeint(mixed_culture_ode, u0, t, args=(p, Competitor_Type),
                         full_output=True, rtol=1e-6, atol=1e-6, tfirst=True)
    return solution


def calculate_doseResponse(u0, t, p, numPoints, Competitor_Type='Resistant'):
    """
        Calculate the change in cell density with different drug
        concentrations. Returns the dose response for each type of 
        microbe, and the range of drug concentrations used.
        
        The method will solve the model `f' using different drug concentrations
        and then retrieve the last density (at t[-1]). Then plot these final
        densities for each drug concentration used.
        
        The number of drug concentrations is determined by `numPoints' and the
        highest concentration by u0[-1].
        
        Assumptions: two types of microbe (S and R).
    """
    drugSupply = u0[-1]
    concentrations_range = np.linspace(0, drugSupply, num=numPoints)
    drS = list()
    drR = list()
    CiS = list()
    CiR = list()
    DiS = list()
    DiR = list()
    for drugSupply in concentrations_range:
        u0[-1] = drugSupply
        sol = ode_solver(u0, t, p, Competitor_Type)
        drS.append(sol.T[0][-1])
        drR.append(sol.T[1][-1])
        CiS.append(sol.T[3][-1])
        CiR.append(sol.T[4][-1])
        DiS.append(sol.T[5][-1])
        DiR.append(sol.T[6][-1])
    return drS, drR, DiS, DiR, CiS, CiR, concentrations_range


def collate_results(solutions, A_range, Inhibition, Competitor_Type=None,
                    Reference='Monoculture', Reference_MIC=None):
    """
        Extract MICs and DICs from dose-response data
    """
    
    S = list()
    R = list()
    MIC_S = list()
    MIC_R = list()
    Carbon_in_S = list()
    Carbon_in_R = list()
    CiS = list()
    CiR = list()
    Drug_in_S = list()
    Drug_in_R = list()
    DiS = list()
    DiR = list()
    for sol in solutions:
        S.append(sol[0].get()[0])  # .get(0) = S, 1 R, 2 A_range
        R.append(sol[0].get()[1])  # .get(0) = S, 1 R, 2 A_range
        Drug_in_S.append(sol[0].get()[2])
        Drug_in_R.append(sol[0].get()[3])
        Carbon_in_S.append(sol[0].get()[4])
        Carbon_in_R.append(sol[0].get()[5])
    S = np.array(S)
    R = np.array(R)
    Drug_in_S = np.array(Drug_in_S)
    Drug_in_R = np.array(Drug_in_R)
    Carbon_in_S = np.array(Carbon_in_S)
    Carbon_in_R = np.array(Carbon_in_R)
    
    if Reference == None:
        # Calculate IC and DiC for the S-type in monoculture conditions.
        for drS, DrugContent, CarbonContent in zip(S, Drug_in_S, Carbon_in_S):
            IC_Calc = interp1d(drS, A_range, bounds_error=False)
            MIC_S.append(IC_Calc(drS[0] * Inhibition).tolist())
            # Drug per cell
            DiC_Calc = interp1d(A_range, DrugContent/drS, bounds_error=False)
            DiS.append(DiC_Calc(IC_Calc(drS[0] * Inhibition)).tolist())
            # Carbon per cell @ MIC
            CiS.append(_calculate_carbon(CarbonContent, A_range, drS,
                                         IC_Calc(drS[0] * Inhibition).tolist()))
    else:
        # Calculate IC for the S-type, and DiC using the MIC in monoculture
        for MIC_ID, zData in enumerate(zip(S, Drug_in_S, Carbon_in_S)):
            drS, DrugContent, CarbonContent = zData
            IC_Calc = interp1d(drS, A_range, bounds_error=False)
            MIC_S.append(IC_Calc(drS[0] * Inhibition).tolist())
            # Drug per cell
            DiC_Calc = interp1d(A_range, DrugContent/drS, bounds_error=False)
            DiS.append(DiC_Calc(Reference_MIC[MIC_ID]).tolist())
            # Carbon per cell @ MIC in Monoculture
            CiS.append(_calculate_carbon(CarbonContent, A_range, drS,
                                         IC_Calc(drS[0] * Inhibition).tolist()))
    
    # Calculate IC for the R-type (only if sensitive)
    if Competitor_Type is not None and Competitor_Type == 'Sensitive':
        for drR, DrugContent, CarbonContent in zip(R, Drug_in_R, Carbon_in_R):
            IC_Calc = interp1d(drR, A_range, bounds_error=False)
            MIC_R.append(IC_Calc(drR[0] * Inhibition).tolist())
            # Drug per cell
            DiC_Calc = interp1d(A_range, DrugContent/drR, bounds_error=False)
            DiR.append(DiC_Calc(IC_Calc(drR[0] * Inhibition)).tolist())
            # Carbon per cell @ MIC in Monoculture
            CiR.append(_calculate_carbon(CarbonContent, A_range, drR,
                                         IC_Calc(drR[0] * Inhibition).tolist()))
    
    return S, MIC_S, R, MIC_R, DiS, DiR, CiS, CiR
