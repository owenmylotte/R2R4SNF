import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import os

# Set up LaTeX fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

def read_horizontal_data(filename):
    """
    Reads horizontally oriented data where each row represents a different measurement
    and columns are different parameters.
    """
    with open(filename, 'r') as f:
        header = f.readline().strip().split(',')
        data_rows = []
        for line in f:
            try:
                row = [float(val) for val in line.strip().split(',')]
                data_rows.append(row)
            except ValueError:
                continue
        data = np.array(data_rows)
        days = data[:, 0]
    return header, days, data


def read_vertical_data(filename):
    """
    Reads vertically oriented data where first column contains labels
    and subsequent columns contain time series data.
    """
    with open(filename, 'r') as f:
        time_points = f.readline().strip().split(',')
        isotopes = []
        data_rows = []
        for line in f:
            values = line.strip().split(',')
            isotope = values[0]
            try:
                row_data = [float(val) if val.strip() else np.nan for val in values[1:]]
                isotopes.append(isotope)
                data_rows.append(row_data)
            except ValueError:
                continue
        activities = np.array(data_rows)
    return time_points, isotopes, activities


def get_top_contributors(days, activities, isotopes, start_idx, end_idx, n_top=10):
    """
    Calculate time-integrated activities using scipy integrator and return top contributors
    """
    integrated_activities = np.zeros(len(isotopes))

    for i, isotope_activities in enumerate(activities):
        valid_activities = isotope_activities[start_idx:end_idx]
        valid_days = days[start_idx:end_idx]

        if not np.any(np.isnan(valid_activities)):
            integrated_activities[i] = integrate.simps(valid_activities, valid_days)

    top_indices = np.argsort(integrated_activities)[-n_top:]
    return top_indices[::-1]

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap, ScalarMappable
from matplotlib.colors import to_rgba

def create_activity_plots(days, activities, isotopes, split_timestep, output_prefix, reactor_name,
                          early_ylim=None, late_ylim=None, colormap_name='gist_rainbow'):
    """
    Creates early and late time plots for a given reactor's data with highly distinguishable isotope colors.
    """
    # Use a perceptually distinct and broad colormap
    cmap = get_cmap(colormap_name)
    num_isotopes = len(isotopes)
    color_map = {isotope: cmap(i / max(num_isotopes - 1, 1)) for i, isotope in enumerate(isotopes)}

    def style_plot(fig, ax):
        ax.set_facecolor('#FFFFFF')  # Background color
        fig.patch.set_facecolor('#FFFFFF')
        ax.grid(True, which="both", ls="-", alpha=0.2, color='#333333')
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color('#CCCCCC')
        ax.tick_params(colors='#333333', labelsize=12)
        ax.xaxis.label.set_color('#333333')
        ax.yaxis.label.set_color('#333333')
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)

    # Find index corresponding to split_timestep
    split_index = np.searchsorted(days, split_timestep)

    # Get top contributors for each period
    top_early = get_top_contributors(days, activities, isotopes, 0, split_index)
    top_late = get_top_contributors(days, activities, isotopes, split_index, len(days))

    # Calculate total activity at each time point
    total_activity = np.nansum(activities, axis=0)

    # Create early-time plot
    fig, ax = plt.subplots(figsize=(12, 7))

    for idx in top_early:
        isotope_activities = activities[idx, :split_index]
        valid_indices = ~np.isnan(isotope_activities)
        if np.any(valid_indices):
            ax.plot(days[:split_index], isotope_activities, label=isotopes[idx],
                    linewidth=1.5, color=color_map[isotopes[idx]])

    # Plot total activity with black dashed line
    ax.plot(days[:split_index], total_activity[:split_index], label='Total Activity',
            linewidth=2.5, color='black', linestyle='--')

    ax.set_xlabel('Time (days)', fontsize=14)
    ax.set_ylabel('Activity', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if early_ylim:
        ax.set_ylim(early_ylim)

    style_plot(fig, ax)

    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,
                    fontsize=12, framealpha=1)
    leg.get_frame().set_facecolor('#FFFFFF')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_early.png', bbox_inches='tight', dpi=500,
                facecolor='#FFFFFF', edgecolor='none')
    plt.close()

    # Create late-time plot
    fig, ax = plt.subplots(figsize=(12, 7))

    for idx in top_late:
        isotope_activities = activities[idx, split_index:]
        valid_indices = ~np.isnan(isotope_activities)
        if np.any(valid_indices):
            ax.plot(days[split_index:], isotope_activities, label=isotopes[idx],
                    linewidth=1.5, color=color_map[isotopes[idx]])

    ax.plot(days[split_index:], total_activity[split_index:], label='Total Activity',
            linewidth=2.5, color='black', linestyle='--')

    ax.set_xlabel('Time (days)', fontsize=14)
    ax.set_ylabel('Activity', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if late_ylim:
        ax.set_ylim(late_ylim)

    style_plot(fig, ax)

    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,
                    fontsize=12, framealpha=1)
    leg.get_frame().set_facecolor('#FFFFFF')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_late.png', bbox_inches='tight', dpi=500,
                facecolor='#FFFFFF', edgecolor='none')
    plt.close()


def convert_half_life_to_decay_constant(half_life, unit='years'):
    """
    Convert half-life to decay constant in per day

    Parameters:
    -----------
    half_life : float
        The half-life value
    unit : str
        The unit of the half-life. Can be 'years', 'days', 'hours', or 'minutes'

    Returns:
    --------
    float : Decay constant in per day
    """
    ln2 = 0.693147180559945

    # Convert all times to days first
    if unit == 'years':
        days = half_life * 365.25
    elif unit == 'days':
        days = half_life
    elif unit == 'hours':
        days = half_life / 24.0
    elif unit == 'minutes':
        days = half_life / (24.0 * 60.0)
    else:
        raise ValueError(f"Unsupported time unit: {unit}")

    return ln2 / days


def get_fissile_data():
    """
    Returns dictionary of fissile isotopes and their decay constants (in per day)
    Half-lives are stored with their units and converted to decay constants
    """
    # Dictionary of isotope: (half_life, unit)
    half_lives = {
        'U232': (68.9, 'years'),
        'U233': (1.59e5, 'years'),
        'U234': (2.46e5, 'years'),
        'U235': (7.04e8, 'years'),
        'U236': (2.34e7, 'years'),
        'U237': (6.75, 'days'),
        'U238': (4.47e9, 'years'),
        'U239': (23.45, 'minutes'),
        'U240': (14.1, 'hours'),
        'U241': (16.0, 'minutes'),
        'Np235': (396.1, 'days'),
        'Np236': (1.55e5, 'years'),
        'Np237': (2.14e6, 'years'),
        'Np238': (2.10, 'days'),
        'Np239': (2.36, 'days'),
        'Pu236': (2.858, 'years'),
        'Pu237': (45.64, 'days'),
        'Pu238': (87.7, 'years'),
        'Pu239': (2.41e4, 'years'),
        'Pu240': (6562.2, 'years'),
        'Pu241': (14.329, 'years'),
        'Pu242': (3.75e5, 'years'),
        'Pu243': (4.955, 'hours'),
        'Am241': (432.6, 'years'),
        'Am242': (16.02, 'hours'),
        'Am243': (7345, 'years'),
        'Am244': (10.2, 'hours'),
        'Am245': (2.05, 'hours'),
        'Cm241': (32.8, 'days'),
        'Cm242': (162.8, 'days'),
        'Cm243': (29.18, 'years'),
        'Cm244': (18.112, 'years'),
        'Cm245': (8245, 'years'),
        'Cm246': (4757, 'years'),
        'Cm247': (1.56e7, 'years'),
        'Cm248': (3.48e5, 'years'),
        'Cm249': (64.15, 'minutes'),
        'Cm250': (8300, 'years')
    }

    # Convert all half-lives to decay constants
    decay_constants = {
        isotope: convert_half_life_to_decay_constant(half_life, unit)
        for isotope, (half_life, unit) in half_lives.items()
    }

    return decay_constants


def create_fissile_inventory_plot(days, inventories, isotopes, split_timestep, output_prefix, reactor_name,
                                  colormap_name='Spectral', early_ylim=None):
    """
    Creates a plot of fissile isotope inventories with perceptually distinct colors,
    split into early and late timesteps.

    Parameters:
        days (array-like): Time points (in days).
        inventories (2D array): Inventory data for each isotope (rows correspond to isotopes, columns to time points).
        isotopes (list of str): Names of the isotopes.
        split_timestep (int): Index at which to split early and late timesteps.
        output_prefix (str): Prefix for the output file name.
        reactor_name (str): Name of the reactor for the plot title.
        colormap_name (str): Name of the matplotlib colormap to use for isotope colors.
        early_ylim (tuple): Optional y-axis limits for the early timestep plot.
    """
    # Use a perceptually distinct colormap
    cmap = get_cmap(colormap_name)
    num_isotopes = len(isotopes)
    color_map = {isotope: cmap(i / max(num_isotopes - 1, 1)) for i, isotope in enumerate(isotopes)}

    # Styling function for plot aesthetics
    def style_plot(fig, ax):
        ax.set_facecolor('#FFFFFF')  # Background color
        fig.patch.set_facecolor('#FFFFFF')
        ax.grid(True, which="both", ls="-", alpha=0.2, color='#333333')
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color('#CCCCCC')
        ax.tick_params(colors='#333333', labelsize=12)
        ax.xaxis.label.set_color('#333333')
        ax.yaxis.label.set_color('#333333')
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)

    # Split the data into early and late timesteps
    early_days = days[:split_timestep]
    late_days = days[split_timestep:]
    early_inventories = inventories[:, :split_timestep]
    late_inventories = inventories[:, split_timestep:]

    # Plot early timestep data
    fig, ax = plt.subplots(figsize=(12, 7))
    for idx, isotope in enumerate(isotopes):
        isotope_inventory = early_inventories[idx, :]
        valid_indices = ~np.isnan(isotope_inventory)
        if np.any(valid_indices):
            ax.plot(early_days, isotope_inventory, label=isotope,
                    linewidth=1.5, color=color_map[isotope])

    total_early_inventory = np.nansum(early_inventories, axis=0)
    ax.plot(early_days, total_early_inventory, label='Total Inventory',
            linewidth=2.5, color='black', linestyle='--')

    ax.set_xlabel('Time (days)', fontsize=14)
    ax.set_ylabel('Inventory (kg)', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if early_ylim:
        ax.set_ylim(early_ylim)
    ax.set_title(f'{reactor_name} Fissile Inventory (Early Time)', fontsize=16, color='#333333')
    style_plot(fig, ax)
    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=12, framealpha=1)
    leg.get_frame().set_facecolor('#FFFFFF')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_early.png', bbox_inches='tight', dpi=500, facecolor='#FFFFFF', edgecolor='none')
    plt.close()

    # Plot late timestep data
    fig, ax = plt.subplots(figsize=(12, 7))
    for idx, isotope in enumerate(isotopes):
        isotope_inventory = late_inventories[idx, :]
        valid_indices = ~np.isnan(isotope_inventory)
        if np.any(valid_indices):
            ax.plot(late_days, isotope_inventory, label=isotope,
                    linewidth=1.5, color=color_map[isotope])

    total_late_inventory = np.nansum(late_inventories, axis=0)
    ax.plot(late_days, total_late_inventory, label='Total Inventory',
            linewidth=2.5, color='black', linestyle='--')

    ax.set_xlabel('Time (days)', fontsize=14)
    ax.set_ylabel('Inventory (kg)', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title(f'{reactor_name} Fissile Inventory (Late Time)', fontsize=16, color='#333333')
    style_plot(fig, ax)
    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=12, framealpha=1)
    leg.get_frame().set_facecolor('#FFFFFF')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_late.png', bbox_inches='tight', dpi=500, facecolor='#FFFFFF', edgecolor='none')
    plt.close()


# Define paths for both reactor types
htgr_base_path = "../../Results/HTGR_FCM/reactor_simulation/"
pwr_base_path = "/Users/owen/PycharmProjects/R2R4SNF/Results/Ref_PWR/reactor_simulation/"
output_path = "../output/"

# Create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

# Process HTGR FCM data
htgr_days_header, htgr_days, htgr_days_data = read_horizontal_data(htgr_base_path + 'summary/out_HTGR_FCM_time_dep.csv')
htgr_time_points, htgr_isotopes, htgr_activities = read_vertical_data(
    htgr_base_path + 'SNF_by_nuclide/radioactivity/HTGR_FCM_activity.csv')

# Process PWR data
pwr_days_header, pwr_days, pwr_days_data = read_horizontal_data(pwr_base_path + 'summary/out_Ref_PWR_time_dep.csv')
pwr_time_points, pwr_isotopes, pwr_activities = read_vertical_data(
    pwr_base_path + 'SNF_by_nuclide/radioactivity/Ref_PWR_activity.csv')

# Create plots with different thresholds and y-axis limits for each reactor type
htgr_split_timestep = 3850
pwr_split_timestep = 1250

# Example y-axis limits (adjust these as needed)
htgr_early_ylim = (1e0, 1e13)
htgr_late_ylim = (1e0, 1e12)
pwr_early_ylim = (1e7, 1e13)
pwr_late_ylim = (1e0, 1e12)

fissile_ylim = (1e-3, 1e0)

create_activity_plots(htgr_days, htgr_activities, htgr_isotopes,
                      htgr_split_timestep, output_path + 'htgr_activity_plot', 'HTGR FCM',
                      early_ylim=htgr_early_ylim, late_ylim=htgr_late_ylim)

create_activity_plots(pwr_days, pwr_activities, pwr_isotopes,
                      pwr_split_timestep, output_path + 'pwr_activity_plot', 'Reference PWR',
                      early_ylim=pwr_early_ylim, late_ylim=pwr_late_ylim)

create_fissile_inventory_plot(htgr_days, htgr_activities, htgr_isotopes, htgr_split_timestep,
                               output_path + 'htgr_fissile_plot', 'HTGR FCM',
                               early_ylim=fissile_ylim)

create_fissile_inventory_plot(pwr_days, pwr_activities, pwr_isotopes, pwr_split_timestep,
                               output_path + 'pwr_fissile_plot', 'Reference PWR',
                               early_ylim=fissile_ylim)
