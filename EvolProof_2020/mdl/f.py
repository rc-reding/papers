def mixed_culture_ode(t, y, p, Competitor_Type='Resistant'):
    """
        Model in Reding (2020)
    """
    # Extract variables
    S, R, C_ext, C_Sint, C_Rint, A_Sint, A_Rint, A_ext = y
    # Extract parameters
    Vmax_S, Vmax_R, Km_S, Km_R,\
        yield_S, yield_R, d, Ki, phi_S, phi_R, Epsilon = p
    
    ### Helper functions ###
    Inhibition_S = 1 / ((1 + A_Sint / Ki)**2)
    Inhibition_R = 1 / ((1 + A_Rint / Ki)**2)
    Uptake_S = Vmax_S * C_ext / (Km_S + C_ext) * Inhibition_S
    
    if Competitor_Type == 'Resistant':
        Uptake_R = Vmax_R * C_ext / (Km_R + C_ext) * (1 - Epsilon)
    elif Competitor_Type == 'Sensitive':
        Uptake_R = Vmax_R * C_ext / (Km_R + C_ext) * Inhibition_R
    
    # To _monitor_ carbon within cells (as `Uptake` is reset every iteration)
    dC_Sintdt = Uptake_S
    dC_Rintdt = Uptake_R
    
    ### Core ODEs ###
    dSdt = Uptake_S * yield_S * S
    dRdt = Uptake_R * yield_R * R
    dC_extdt = - (Uptake_S + Uptake_R) * (S + R)
    
    dA_Sintdt = -d * A_Sint + phi_S * S * (A_ext - A_Sint)
    dA_Rintdt = -d * A_Rint + phi_R * R * (A_ext - A_Rint)
    dA_extdt = -d * A_ext - phi_S * (A_ext - A_Sint) * S - phi_R * (A_ext - A_Rint) * R
    
    return [dSdt, dRdt, dC_extdt, dC_Sintdt, dC_Rintdt,
            dA_Sintdt, dA_Rintdt, dA_extdt]

def set_ode_params(Vmax_S=1.25, Vmax_R=1.25, Km_S=0.5, Km_R=0.5, yield_S=0.65,
                   yield_R=0.65, d=0.0001, Ki=0.1, phi_S=10.0, phi_R=10.0,
                   Epsilon=0.0):
    """
        Parameters:
        
        - Vmax: Maximal carbon uptake rate.
        - Km: Affinity for carbon source, a.k.a. half-saturation parameter.
        - yield: Biomass yield as given by K/c (Monod 1947), where "K" is the
            population size at the equilibrium and "c" the amount of carbon
            supplied.
        - Epsilon: Costs of resistance.
    """
    return list([Vmax_S, Vmax_R, Km_S, Km_R, yield_S, yield_R,
                  d, Ki, phi_S, phi_R, Epsilon])
