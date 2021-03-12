#!/usr/bin/env python3

from matplotlib import pyplot as plt
import multiprocessing as mp
import numpy as np

from mdl import set_ode_params
from helper_func import *


# Set values for initial conditions vector u0:
A_Sint = 0.0
A_Rint = 0.0
A_init = 2.0

C_Sint = 0.0
C_Rint = 0.0
C_init = 2.0

S_init = 0.001
R_init = 0.001

u0 = [S_init, R_init, C_init, C_Sint, C_Rint, A_Sint, A_Rint, A_init]

# Generate parameters vector:
t = np.linspace(0, 24, num=100)
p = set_ode_params()

numPoints = 150
InhibitionFactor = 0.9
Inhibition = 1 - InhibitionFactor

# Sensitivity of competing species: 'Resistant' or 'Sensitive'.
R_Type = 'Resistant'

# Growth conditions: 'Monoculture' or 'Competition' (default).
R_Condition = 'Competition'

# Generate reference:
A_range = np.linspace(0, A_init, num=numPoints)


### MULTIPROCESSING STEP A: Initialise workers. ###
CPU_N = mp.cpu_count()  # For multiprocessing.
pool = mp.Pool(processes=CPU_N)

u0[1] = 0.0  # Remove competitor.
# Generate reference doseResponse profiles varying Vmax_S:
_, Ref_VmaxS_MIC_S, _, _, _,\
    Ref_VmaxS_DiS, _, Ref_VmaxS_CiS, _ = screen_VmaxS(pool, A_range, u0, t, p,
                                                      Inhibition, numPoints)
u0[1] = R_init  # Introduce competitor.
VmaxS_HeatMap, VmaxS_MIC_S, VmaxR_HeatMap,\
    VmaxR_MIC_R, VmaxS_range, VmaxS_DiS,\
    VmaxS_DiR, VmaxS_CiS, VmaxS_CiR = screen_VmaxS(pool, A_range, u0, t, p,
                                                  Inhibition, numPoints,
                                                  R_Type, R_Condition,
                                                  Ref_VmaxS_MIC_S)

u0[1] = 0.0  # Remove competitor.
# Generate reference doseResponse profiles varying Km_S:
_, Ref_KmS_MIC_S, _, _, _,\
    Ref_KmS_DiS, _, Ref_KmS_CiS, _ = screen_KmS(pool, A_range, u0, t, p,
                                                Inhibition, numPoints)
u0[1] = R_init  # Introduce competitor.
KmS_HeatMap, KmS_MIC_S, KmR_HeatMap, KmR_MIC_R,\
    KmS_range, KmS_DiS, KmS_DiR, KmS_CiS,\
    KmS_CiR = screen_KmS(pool, A_range, u0, t, p, Inhibition, numPoints,
                         R_Type, R_Condition, Ref_KmS_MIC_S)

u0[1] = 0.0  # Remove competitor.
# Generate reference doseResponse profiles varying Yield_S:
_, Ref_YieldS_MIC_S, _, _, _,\
    Ref_YieldS_DiS, _, Ref_YieldS_CiS, _ = screen_YieldS(pool, A_range, u0, t,
                                                         p, Inhibition,
                                                         numPoints)
u0[1] = R_init  # Introduce competitor.
YieldS_HeatMap, YieldS_MIC_S, YieldR_HeatMap,\
    YieldR_MIC_R, YieldS_range, YieldS_DiS,\
    YieldS_DiR, YieldS_CiS, YieldS_CiR = screen_YieldS(pool, A_range, u0, t,
                                                       p, Inhibition,
                                                       numPoints, R_Type,
                                                       R_Condition,
                                                       Ref_YieldS_MIC_S)
u0[1] = 0.0  # Remove competitor.

### MULTIPROCESSING STEP B: Close workers. ###
pool.close()
pool.join()


""" PLOTTING SECTION """

# HEATMAPS
plotHeatMap(A_range, VmaxS_range, VmaxS_HeatMap, VmaxS_MIC_S, "Vmax_S",
            LabelY=r"$\bar{\mathrm{\mu}}_{1}$ (mg/OD/h)",
            Reference_MIC_S=Ref_VmaxS_MIC_S, P_Competitor=p[1],
            Competitor_Type=R_Type)
plotHeatMap(A_range, KmS_range, KmS_HeatMap, KmS_MIC_S, "Km_S",
            LabelY="k$_1$ (mg/mL)", Reference_MIC_S=Ref_KmS_MIC_S,
            P_Competitor=p[3], Competitor_Type=R_Type)
plotHeatMap(A_range, YieldS_range, YieldS_HeatMap, YieldS_MIC_S, "Yield_S",
            LabelY="y$_1$ (OD/mg)", Reference_MIC_S=Ref_YieldS_MIC_S,
            P_Competitor=p[5], Competitor_Type=R_Type)

# DIFFERENCES IN MIC
plot_IC_differences(VmaxS_range, Ref_VmaxS_MIC_S, VmaxS_MIC_S, "Vmax_S",
                    LabelX=r"$\bar{\mathrm{\mu}}_{1}$ (mg/OD/h)",
                    Competitor_Type=R_Type)
plot_IC_differences(KmS_range, Ref_KmS_MIC_S, KmS_MIC_S, "Km_S",
                    LabelX="k$_1$ (mg/mL)", Competitor_Type=R_Type)
plot_IC_differences(YieldS_range, Ref_YieldS_MIC_S, YieldS_MIC_S, "Yield_S",
                    LabelX="y$_1$ (OD/mg)", Competitor_Type=R_Type)

# Plot drug per cell of focal species S, for each trait changed
plot_p_vs_DiffContent(Ref_VmaxS_MIC_S, list([Ref_VmaxS_DiS, VmaxS_DiS]),
                    "Vmax_S", R_Type, Content="Drug",
                    LabelX=r"IC$_{90}^*$ (\textmu g/mL), $\bar{\mu}_{1}$")
plot_p_vs_DiffContent(Ref_KmS_MIC_S, list([Ref_KmS_DiS, KmS_DiS]), "Km_S",
                    R_Type, Content="Drug",
                    LabelX=r"IC$_{90}^*$ (\textmu g/mL), k$_1$")
plot_p_vs_DiffContent(Ref_YieldS_MIC_S, list([Ref_YieldS_DiS, YieldS_DiS]),
                    "Yield_S", R_Type, Content="Drug",
                    LabelX=r"IC$_{90}^*$ (\textmu g/mL), y$_1$")

# Plot explanatory 3D heat map.
if R_Type == 'Resistant':
    plot3D_to_2D(A_range, YieldS_range, VmaxS_HeatMap, VmaxS_MIC_S, R_Type)
