import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker as tkr
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch


class Arrow3D(FancyArrowPatch):
    """
        Plotting arrow in 3D space. Snippet from:
        https://stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d

        Args:
        xs: vector containing initial and final position of the arrow in the X axis.
        ys: vector containing initial and final position of the arrow in the Y axis.
        zs: vector containing initial and final position of the arrow in the Z axis.
        *args: check matplotlib documentation.
        **kwargs: check matplotlib documentation.
    """
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


def plotHeatMap(A_range, P_range, P_HeatMap, MIC_S, P_Label,
                LabelX=r"Antibiotic (\textmu g/mL)", LabelY=None,
                Reference_MIC_S=None, P_Competitor=None,
                Competitor_Type='Resistant'):
    """
        Represent cell density at different drug concentrations AND with
        different parameter values.

        Args:
        A_range: Range of antimicrobial (drug) used. Array-type.
        P_range: Range of values for the parameter `p' tweaked.
        P_HeatMap: MxN matrix containing cell density.
        MIC_S: Minimum Inhibitory Inhibition of the most sensitive type of
        microbe measured in mixed culture conditions.
        P_Label: Name of the parameter tested.
        Label_X: Label used for the X axis of the plot.
        Label_Y: Label used for the Y axis of the plot.
        Reference_MIC_S: Minimum Inhibitory Inhibition of the most sensitive
        type of microbe measured in pure culture conditions.
        P_competitor: Parameter value for the competing species, used as a
        reference.
    """
    fig = plt.figure(figsize=(5, 5))
    ax = plt.gca()  # Create axis (only if I need to manipulate it, and I do).
    ax.tick_params(axis='both', direction='out', labelsize=20)

    plt.pcolormesh(A_range, P_range, P_HeatMap, cmap='bone_r', linewidth=0,
                   shading='auto')
    plt.plot(MIC_S, P_range, linewidth=3, color='black', alpha=0.5)

    if P_Competitor:
        plt.plot(A_range[-1] * -0.005, P_Competitor, marker=5, markersize=6,
                 color='black', alpha=0.05, clip_on=False)
        plt.plot([A_range[0], A_range[-1]], [P_Competitor, P_Competitor],
                 linewidth=1, linestyle='dashed', dashes=(5, 5),
                 color='black')
        if P_Label == 'Km_S' and Competitor_Type == 'Sensitive':
            plt.text(1.5, 0.4, 'Parameter value\nfor species S$_2$',
                     fontsize=18, verticalalignment='center',
                     horizontalalignment='center')

    plt.xlabel(LabelX, fontsize=24)
    plt.ylabel(LabelY, fontsize=24)

    plt.xlim([A_range[0], A_range[-1]])

    fx = lambda x, pos: str(x).rstrip('0').rstrip('.')  # int shown as int when float are in axis.
    fy = lambda y, pos: str(round(y, ndigits=4)).rstrip('0').rstrip('.')  # int shown as int when float are in axis.
    ax.xaxis.set_major_formatter(tkr.FuncFormatter(fx))
    ax.yaxis.set_major_formatter(tkr.FuncFormatter(fy))

    if Reference_MIC_S:
        plt.plot(Reference_MIC_S, P_range, linewidth=3, color='darkgray',
                 alpha=0.5)
        plt.savefig('./img/' + P_Label + '_Competition_' + Competitor_Type +
                    '.eps', bbox_inches='tight', transparent=True)
        plt.close()
    else:
        plt.savefig('./img/' + P_Label + '.eps', bbox_inches='tight',
                    transparent=True)
        plt.close()


def plot3D_to_2D(A_range, P_range, P_HeatMap, P_MIC,
                 Competitor_Type='Resistant'):
    """ Surface plot overlaying a mesh plot to summarise the meaning of the
        heat maps used.
    """
    A_grid, P_grid = np.meshgrid(A_range, P_range)
    Z_range = np.ones_like(A_range) * 0.1  # Height for MIC line
    fig = plt.figure(figsize=(10, 8))
    ax = fig.gca(projection='3d')  # Create axis (only if I need to manipulate it, and I do).
    ax.tick_params(axis='both', direction='out', labelsize=11)
    # Plot surface first:
    surf_plot = ax.plot_surface(A_grid, P_grid, P_HeatMap, linestyle='None',
                                shade=False, edgecolor='none', cmap='bone_r',
                                ccount=A_range.shape[0]/2,
                                rcount=P_range.shape[0])
    ax.plot_wireframe(A_grid, P_grid, P_HeatMap, linewidth=0.75,
                      ccount=0, rcount=15, zorder=5)
    ax.plot3D(P_MIC, P_range, Z_range, color='black', zorder=4,
              linewidth=2.25)
    # Annotations
    ax.text(1.25, 1.75, 1.2, 'Dose-response with\nparameter value $n$',
            fontsize=12, verticalalignment='center',
            horizontalalignment='center', backgroundcolor='none')
    arr3D = Arrow3D([1, 0.7], [2, 2], [0.95, 0.35], mutation_scale=15,
                    linewidth=1.25, arrowstyle='-|>', color='k')
    ax.add_artist(arr3D)
    ax.text(1.65, 1.75, 0.6, 'IC$_{90}$ with\nparameter value $n$',
            fontsize=12, verticalalignment='center',
            horizontalalignment='center', backgroundcolor='none')
    arr3D = Arrow3D([1.35, 1.175], [2.25, 2.25], [0.25, 0], mutation_scale=15,
                    linewidth=1.25, arrowstyle='-|>', color='k')
    ax.add_artist(arr3D)
    # Axis style
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis._axinfo['grid'].update({'linewidth': 0.25, 'color': 'gainsboro'})
    ax.yaxis._axinfo['grid'].update({'linewidth': 0.25, 'color': 'gainsboro'})
    ax.zaxis._axinfo['grid'].update({'linewidth': 0.25, 'color': 'gainsboro'})
    # Axis labels
    ax.set_ylabel('Parameter Values', fontsize=16, labelpad=10)
    ax.set_xlabel('Antibiotic Concentration', fontsize=16, labelpad=10)
    ax.set_zlabel('Cell Density', fontsize=16, labelpad=10)
    fx = lambda x, pos: str(x).rstrip('0').rstrip('.')  # int shown as int when float are in axis.
    ax.xaxis.set_major_formatter(tkr.FuncFormatter(fx))
    ax.yaxis.set_major_formatter(tkr.FuncFormatter(fx))
    # Colorbar
    cb = fig.colorbar(surf_plot, ax=ax, ticks=[0.1, 0.9],
                      orientation='horizontal', shrink=0.25, pad=0.075)
    cb.ax.set_xticklabels(['Low', 'High'])
    cb.set_label('Cell Density', size=12)
    # Save
    fig.savefig('./img/3D_HeatMap_' + Competitor_Type + '.eps',
                bbox_inches='tight', pad_inches=0.3)
    plt.close()


def plot_IC_differences(P_range, Reference_MIC_S, MIC_S, P_Label, LabelX=None,
                        Competitor_Type='Resistant'):
    """ Plot difference in MICs between two conditions """
    delta_MIC = np.array(MIC_S) - np.array(Reference_MIC_S)

    fig = plt.figure(figsize=(5, 5))
    ax = plt.gca()  # Create axis (only if I need to manipulate it, and I do).
    ax.tick_params(axis='both', which='both', direction='in', labelsize=20)
    # Plot
    plt.plot(P_range, np.zeros_like(P_range), color='black', alpha=0.75,
             linewidth=0.75, linestyle='dashed', dashes=(5, 5))
    plt.plot(P_range, delta_MIC, color='black', linewidth=4)
    # Axis
    f = lambda x, pos: str(round(x, ndigits=4)).rstrip('0').rstrip('.')  # int shown as int when float are in axis.
    ax.yaxis.set_major_formatter(tkr.FuncFormatter(f))
    ax.xaxis.set_major_formatter(tkr.FuncFormatter(f))
    
    ax.set_xlim(P_range.min() * 0.95, P_range.max() * 1.05)
    # Labels
    plt.xlabel(LabelX, fontsize=22)
    plt.ylabel("Difference in IC$_{90}$\n" + r"S$_1$ (\textmu g/mL)",
               fontsize=22)
    # Save
    plt.savefig('./img/MIC_Difference_' + P_Label + '_' + Competitor_Type +
                '.eps', bbox_inches='tight')
    plt.close()


def plot_p_vs_DiffContent(P_range, DiCs, P_Label, Competitor_Type,
                          Content='Drug', LabelX='Drug', Competition=True):
    """
        Plot antimicrobial content per cell (drug inside cells, or DIC)
        against parameter range in monoculture and, optionally, mixed culture
        conditions.
    """
    P_range = np.array(P_range)
    
    if Competition is True:
        DiC_mono = DiCs[0]
        DiC_mixed = DiCs[1]
    else:
        DiC_mono = DiCs
    
    if Content == 'Drug':
        F_Label = "DIC"
        Content_Units = r"\textmu g/mL/cell"
    elif Content == 'Carbon':
        F_Label = "Carbon"
        Content_Units = "mg/mL/cell"
    
    fig = plt.figure(figsize=(5, 5))
    ax = plt.gca()  # Create axis (only if I need to manipulate it, and I do).
    ax.tick_params(axis='both', which='both', direction='in', labelsize=20)
    # Plot
    if Competition is True:
        plt.plot([P_range.min(), P_range.max()], np.zeros_like([P_range.min(),
                 P_range.max()]), color='black', alpha=0.75, linewidth=0.75,
                 linestyle='dashed', dashes=(5, 5))
        plt.plot(P_range, np.array(DiC_mono) - np.array(DiC_mixed),
                 color='black', alpha=0.5, linewidth=3)
        # Axes
        plt.ylabel('Difference in relative\n' + Content.lower() +
                   ' content (' + Content_Units + ')', fontsize=24)
        if P_Label == 'Km_S' and Competitor_Type == 'Sensitive':
            # Arrow
            ax.annotate('', xy=(0.625, 0.75), xytext=(0.625,-0.75),
                        arrowprops=dict(arrowstyle='<|-|>'))
            plt.text(0.51, 0.25, 'More drug in pure culture',
                     fontsize=16, verticalalignment='center',
                     horizontalalignment='left')
            plt.text(0.51, -0.25, 'More drug in mixed culture',
                     fontsize=16, verticalalignment='center',
                     horizontalalignment='left')
            # Axes
            ylims = ax.get_ylim()
            ax.set_ylim(-1.2, ylims[1])
    else:
        plt.plot(P_range, DiC_mono, color='darkgray', alpha=0.75, linewidth=3)
        plt.ylabel(Content + ' per S-cell\n(' + Content_Units + ')',
                   fontsize=24)
        # Axes
        ax.grid(True, which='both', alpha=0.25, linestyle='dotted')
    
    # Axes
    plt.xlabel(LabelX, fontsize=22)
    if P_Label == 'Km_S':
        ax.set_xlim(P_range.min() * 0.99975, P_range.max() * 1.00025)
    else:
        ax.xaxis.set_major_locator(tkr.MultipleLocator(0.2))
        ax.set_xlim(P_range.min() * 0.95, P_range.max() * 1.01)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    
    f = lambda x, pos: str(round(x, ndigits=4)).rstrip('0').rstrip('.')  # int shown as int when float are in axis.
    ax.yaxis.set_major_formatter(tkr.FuncFormatter(f))
    ax.xaxis.set_major_formatter(tkr.FuncFormatter(f))
    
    # Save
    plt.savefig('./img/' + P_Label + '_vs_' + F_Label + '_' +
                Competitor_Type + '.eps', bbox_inches='tight')
    plt.close()

