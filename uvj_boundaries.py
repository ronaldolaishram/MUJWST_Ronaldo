"""
UVJ Diagram Boundary Functions
Williams et al. 2009 and Whitaker et al. 2011 Selection Criteria
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path


def UVJ_separation(UV, VJ, redshift=2.1):
    """
    Whitaker+ 11 colour-cut criteria with redshift dependence

    Returns:
    --------
    int : 1 for quiescent, 0 for star-forming, -1 for undefined redshift
    """
    if redshift > 4:
        return -1
    if redshift < 0.5:
        if (UV > 0.88*VJ+0.69 and UV > 1.3 and VJ < 1.6):
            return 1
        else:
            return 0
    if redshift >= 0.5 and redshift < 1.5:
        if (UV > 0.88*VJ+0.59 and UV > 1.3 and VJ < 1.6):
            return 1
        else:
            return 0
    if redshift >= 1.5 and redshift < 2.0:
        if (UV > 0.88*VJ+0.59 and UV > 1.3 and VJ < 1.5):
            return 1
        else:
            return 0
    if redshift >= 2 and redshift < 4:
        if (UV > 0.88*VJ+0.59 and UV > 1.2 and VJ < 1.4):
            return 1
        else:
            return 0


def is_quiescent_dusty(uv, vj, redshift=2.16, criteria='williams'):
    """
    Determine if a pixel is quiescent or star-forming based on UVJ criteria.
    Dusty galaxies are considered star-forming.

    Parameters:
    -----------
    uv : float or array
        U-V color
    vj : float or array
        V-J color
    redshift : float
        Galaxy redshift (affects diagonal boundary)
    criteria : str
        Selection criteria ('williams')

    Returns:
    --------
    int : 2 for quiescent, 0 for star-forming (including dusty)
    """
    if criteria.lower() == 'williams':
        # Apply diagonal boundary check based on redshift
        if 0.0 <= redshift < 0.5:
            above_diagonal = uv > (0.88 * vj + 0.69)
        elif 0.5 <= redshift < 1.0:
            above_diagonal = uv > (0.88 * vj + 0.59)
        elif 1.0 <= redshift < 2.0:
            above_diagonal = uv > (0.88 * vj + 0.49)
        else:  # z >= 2.0 - extrapolate using z=1-2 criteria
            above_diagonal = uv > (0.88 * vj + 0.49)

        # Quiescent: above diagonal AND UV > 1.3 AND VJ < 1.6
        if above_diagonal and uv > 1.3 and vj < 1.6:
            return 2  # Quiescent
        else:
            return 0  # Star-forming (includes dusty)

    else:
        raise ValueError("Only Williams criteria is currently supported")


def is_quiescent_dusty_fixed(uv, vj, redshift=2.16, criteria='williams'):
    """
    Determine if a pixel is quiescent, star-forming, or dusty based on UVJ criteria.
    Fix ensures proper quiescent classification strictly within UVJ boundaries.

    Parameters:
    -----------
    uv : float or array
        U-V color
    vj : float or array
        V-J color
    redshift : float
        Galaxy redshift
    criteria : str
        Selection criteria ('williams')

    Returns:
    --------
    int : 2 for quiescent, 1 for dusty star-forming, 0 for star-forming
    """
    if criteria.lower() == 'williams':
        # Define the polygon for the quiescent region
        quiescent_region = Path([
            (0.0, 1.3), (0.85, 1.3), (1.6, 1.95), (1.6, 2.5), (0.0, 2.5)
        ])

        # Define the polygon for the dusty region
        dusty_region = Path([
            (1.2, 0.0), (1.2, 1.60333), (1.6, 1.95), (1.6, 2.5), (2.5, 2.5), (2.5, 0.0)
        ])

        # Check if the pixel is in the quiescent region
        if quiescent_region.contains_point((vj, uv)):
            return 2  # Quiescent

        # Check if the pixel is in the dusty region
        if dusty_region.contains_point((vj, uv)):
            return 1  # Dusty star-forming

        # Otherwise, classify as star-forming
        return 0
    else:
        raise ValueError("Only Williams criteria is currently supported")


def plot_uvj_boundaries(ax, add_labels=True, add_shading=True, redshift=2.1):
    """
    Plot UVJ selection boundaries using exact coordinates

    Parameters:
    -----------
    ax : matplotlib axis
        Axis to plot on
    add_labels : bool
        Whether to add text labels for regions
    add_shading : bool
        Whether to shade the regions
    redshift : float
        Galaxy redshift (for reference)

    Returns:
    --------
    None (modifies ax in place)
    """

    # Define UVJ regions - adjusted to match classification function
    # Points on diagonal line at z > 2: UV = 0.88 * VJ + 0.49
    diag_vj = np.array([0.92, 1.6])
    diag_uv = 0.88 * diag_vj + 0.49  # Corresponding UV points

    # Quiescent region polygon
    quiescent_region = [
        [-0.4, 1.3],           # Start at left edge
        [-0.4, 2.5],           # Up to top left
        [1.6, 2.5],            # Top right
        [1.6, diag_uv[1]],     # Down to diagonal intersection
        [diag_vj[0], diag_uv[0]]  # Follow diagonal to left
    ]

    # Star-forming region - everything outside quiescent region
    sf_region = [
        [-0.4, 0],                        # Bottom left
        [-0.4, 1.3],                      # Up to UV limit
        [0.92, 1.3],                      # To where diagonal intersects UV = 1.3
        [0.92, 0.88 * 0.92 + 0.49],       # Follow diagonal up to VJ = 0.92
        [2.3, 0.88 * 2.3 + 0.49],         # Extend to right edge along diagonal
        [2.3, 2.5],                       # Up to top right
        [1.6, 2.5],                       # Left to quiescent boundary
        [1.6, 0.88 * 1.6 + 0.49],         # Down diagonal
        [2.3, 0.88 * 2.3 + 0.49],         # Back to diagonal
        [2.3, 0]                          # Down to bottom
    ]

    # Add polygons for the regions
    if add_shading:
        ax.add_patch(patches.Polygon(quiescent_region, closed=True,
                                    fill=True, facecolor=(1, 0, 0, 0.05),
                                    edgecolor='black', linewidth=1))
        ax.add_patch(patches.Polygon(sf_region, closed=True,
                                    fill=True, facecolor=(0, 0, 1, 0.05)))

    # Add diagonal line (optional - commented out in original)
    # vj_line = np.linspace(-0.4, 2.2, 100)
    # uv_line = 0.88 * vj_line + 0.49
    # ax.plot(vj_line, uv_line, '--', color='black', alpha=0.5)

    # Add text annotations if requested
    if add_labels:
        plt.annotate('Quiescent Region', xy=(-0.3, 2.3), color='black', fontsize=12)
        plt.annotate('Star Forming Region', xy=(-0.3, 0.1), color='black', fontsize=12)


def create_uvj_plot(V_J, U_V, classification_mask=None,
                    figsize=(12, 9), xlim=(-0.4, 2.3), ylim=(0, 2.5),
                    redshift=2.1):
    """
    Create a complete UVJ diagram with selection boundaries
    Uses exact implementation from create_galaxy_uvj_map_dusty

    Parameters:
    -----------
    V_J : array
        V-J colors for all galaxies
    U_V : array
        U-V colors for all galaxies
    classification_mask : array of int, optional
        If provided, 2 = quiescent, 0 = star-forming
        If None, will classify using is_quiescent_dusty
    figsize : tuple
        Figure size
    xlim : tuple
        X-axis limits (V-J)
    ylim : tuple
        Y-axis limits (U-V)
    redshift : float
        Galaxy redshift (affects classification)

    Returns:
    --------
    fig, ax : matplotlib figure and axis
    """

    # Classify if not provided
    if classification_mask is None:
        classification_mask = np.array([is_quiescent_dusty(uv, vj, redshift)
                                       for uv, vj in zip(U_V, V_J)])

    quiescent_mask = (classification_mask == 2)
    star_forming_mask = (classification_mask == 0)

    n_quiescent = np.sum(quiescent_mask)
    n_sf = np.sum(star_forming_mask)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot star-forming galaxies
    ax.scatter(V_J[star_forming_mask], U_V[star_forming_mask],
               c='blue', s=5, alpha=0.3, label=f'Star-forming ({n_sf})', rasterized=True)

    # Plot quiescent galaxies
    ax.scatter(V_J[quiescent_mask], U_V[quiescent_mask],
               c='red', s=10, alpha=0.5, label=f'Quiescent ({n_quiescent})', rasterized=True)

    # Add UVJ boundaries
    plot_uvj_boundaries(ax, add_labels=True, add_shading=True, redshift=redshift)

    # Labels and formatting
    ax.set_xlabel('V-J', fontsize=16)
    ax.set_ylabel('U-V', fontsize=16)
    ax.set_title('UVJ Classification: Star-Forming vs Quiescent', fontsize=18, fontweight='bold')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend(loc='upper left', fontsize=12, framealpha=0.9, markerscale=3)
    # ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, ax


# Example usage:
if __name__ == "__main__":
    # Example: generate random galaxy colors
    np.random.seed(42)
    n_galaxies = 1000

    # Simulate blue cloud (star-forming)
    V_J_sf = np.random.normal(0.5, 0.4, 700)
    U_V_sf = np.random.normal(0.8, 0.3, 700)

    # Simulate red sequence (quiescent)
    V_J_q = np.random.normal(0.8, 0.3, 300)
    U_V_q = np.random.normal(1.7, 0.2, 300)

    # Combine
    V_J = np.concatenate([V_J_sf, V_J_q])
    U_V = np.concatenate([U_V_sf, U_V_q])

    # Create plot
    fig, ax = create_uvj_plot(V_J, U_V, redshift=2.1)
    plt.show()

# %%
