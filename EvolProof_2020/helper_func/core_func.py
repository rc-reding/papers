import numpy as np
from .ds_solver import calculate_doseResponse, collate_results


def screen_VmaxS(pool, A_range, u0, t, p, Inhibition, numPoints,
                 Competitor_Type='Resistant', Reference=None,
                 Reference_MIC=None, p_min=0.35, p_max=5.0):
    """
        Screen through Vmax data of the most sensitive type of microbe (S) and
        return cell density (S), MIC, and parameter range used. Note that S has
        TWO, not one, dimensions: cell density with different drug
        concentrations AND with different values for Vmax.

        Arguments:
        pool: A process pool object which controls a pool of worker processes
        to which jobs can be submitted.
        A_range: Range of antimicrobial (drug) used. Array-type.
        Inhibition: Relative density sought after with the antibiotic. `0'
        denotes absence of growth (100% inhibition) whereas `1' is full growth
        (0% inhibition).
        p_min: Lowest Vmax to be tested.
        p_max: Highest Vmax to be tested.

        NOTE: This method cannot be abstracted due to the use of the
        `multiprocessing' library. Small price to pay for the _huge_ gain in
        performance.
    """
    solutions = list()
    Vmax_S, Vmax_R, Km_S, Km_R,\
        yield_S, yield_R, deg, Ki, phi_S, phi_R, Epsilon = p
    par_range = np.linspace(p[0] * p_min, p[0] * p_max, num=numPoints)
    for par_value in par_range:
        solutions.append([pool.apply_async(calculate_doseResponse, args=(u0,
                                     t, [par_value, Vmax_R, Km_S, Km_R,
                                     yield_S, yield_R, deg, Ki, phi_S, phi_R,
                                     Epsilon], numPoints, Competitor_Type))])
    
    S, MIC_S, R, MIC_R,\
        DiS, DiR, CiS, CiR = collate_results(solutions, A_range, Inhibition,
                                             Competitor_Type, Reference,
                                             Reference_MIC)
    return S, MIC_S, R, MIC_R, par_range, DiS, DiR, CiS, CiR


def screen_KmS(pool, A_range, u0, t, p, Inhibition, numPoints,
               Competitor_Type='Resistant', Reference=None,
               Reference_MIC=None, p_min=0.01, p_max=2.5):
    """
        Screen through Km data of the most sensitive type of microbe (S) and
        return cell density (S), MIC, and parameter range used. Note that S has
        TWO, not one, dimensions: cell density with different drug
        concentrations AND with different values for Km.

        Arguments:
        pool: A process pool object which controls a pool of worker processes
        to which jobs can be submitted.
        A_range: Range of antimicrobial (drug) used. Array-type.
        Inhibition: Relative density sought after with the antibiotic. `0'
        denotes absence of growth (100% inhibition) whereas `1' is full growth
        (0% inhibition).
        p_min: Lowest Km to be tested.
        p_max: Highest Km to be tested.

        NOTE: This method cannot be abstracted due to the use of the
        `multiprocessing' library. Small price to pay for the _huge_ gain in
        performance.
    """
    solutions = list()
    Vmax_S, Vmax_R, Km_S, Km_R,\
        yield_S, yield_R, deg, Ki, phi_S, phi_R, Epsilon = p
    par_range = np.linspace(p[2] * p_min, p[2] * p_max, num=numPoints)
    for par_value in par_range:
        solutions.append([pool.apply_async(calculate_doseResponse, args=(u0,
                                     t, [Vmax_S, Vmax_R, par_value, Km_R,
                                     yield_S, yield_R, deg, Ki, phi_S, phi_R,
                                     Epsilon], numPoints, Competitor_Type))])
    
    S, MIC_S, R, MIC_R,\
        DiS, DiR, CiS, CiR = collate_results(solutions, A_range, Inhibition,
                                             Competitor_Type, Reference,
                                             Reference_MIC)
    return S, MIC_S, R, MIC_R, par_range, DiS, DiR, CiS, CiR


def screen_YieldS(pool, A_range, u0, t, p, Inhibition, numPoints,
                  Competitor_Type='Resistant', Reference=None,
                  Reference_MIC=None, p_min=0.25, p_max=3.0):
    """
        Screen through Yield data of the most sensitive type of microbe (S) and
        return cell density (S), MIC, and parameter range used. Note that S has
        TWO, not one, dimensions: cell density with different drug
        concentrations AND with different values for Yield.

        Arguments:
        pool: A process pool object which controls a pool of worker processes
        to which jobs can be submitted.
        A_range: Range of antimicrobial (drug) used. Array-type.
        Inhibition: Relative density sought after with the antibiotic. `0'
        denotes absence of growth (100% inhibition) whereas `1' is full growth
        (0% inhibition).
        p_min: Lowest Yield to be tested.
        p_max: Highest Yield to be tested.

        NOTE: This method cannot be abstracted due to the use of the
        `multiprocessing' library. Small price to pay for the _huge_ gain in
        performance.
    """
    solutions = list()
    Vmax_S, Vmax_R, Km_S, Km_R,\
        yield_S, yield_R, deg, Ki, phi_S, phi_R, Epsilon = p
    par_range = np.linspace(p[4] * p_min, p[4] * p_max, num=numPoints)
    for par_value in par_range:
        solutions.append([pool.apply_async(calculate_doseResponse, args=(u0,
                                         t, [Vmax_S, Vmax_R, Km_S, Km_R,
                                         par_value, yield_R, deg, Ki, phi_S,
                                         phi_R, Epsilon], numPoints,
                                         Competitor_Type))])
    
    S, MIC_S, R, MIC_R,\
        DiS, DiR, CiS, CiR = collate_results(solutions, A_range, Inhibition,
                                             Competitor_Type, Reference,
                                             Reference_MIC)
    return S, MIC_S, R, MIC_R, par_range, DiS, DiR, CiS, CiR
