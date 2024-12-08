import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import os


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


def create_activity_plots(days, activities, isotopes, split_timestep, output_prefix, reactor_name,
                          early_ylim=None, late_ylim=None):
    """
    Creates early and late time plots for a given reactor's data

    Parameters:
    -----------
    early_ylim : tuple, optional
        (ymin, ymax) for early-time plot
    late_ylim : tuple, optional
        (ymin, ymax) for late-time plot
    """
    # Find index corresponding to split_timestep
    split_index = np.searchsorted(days, split_timestep)

    # Get top contributors for each period
    top_early = get_top_contributors(days, activities, isotopes, 0, split_index)
    top_late = get_top_contributors(days, activities, isotopes, split_index, len(days))

    # Calculate total activity at each time point
    total_activity = np.nansum(activities, axis=0)

    # Create early-time plot
    plt.figure(figsize=(12, 8))
    for idx in top_early:
        isotope_activities = activities[idx, :split_index]
        valid_indices = ~np.isnan(isotope_activities)
        if np.any(valid_indices):
            plt.plot(days[:split_index], isotope_activities, label=isotopes[idx], linewidth=1)

    plt.plot(days[:split_index], total_activity[:split_index],
             label='Total Activity', linewidth=2, color='black', linestyle='--')

    plt.xlabel('Time (days)')
    plt.ylabel('Activity')
    plt.title(f'{reactor_name} In Core Activity\nTop 10 Time Integrated Contributors')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.xscale('log')
    plt.yscale('log')
    if early_ylim:
        plt.ylim(early_ylim)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize='small')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_early.png', bbox_inches='tight', dpi=500)
    plt.close()

    # Create late-time plot
    plt.figure(figsize=(12, 8))
    for idx in top_late:
        isotope_activities = activities[idx, split_index:]
        valid_indices = ~np.isnan(isotope_activities)
        if np.any(valid_indices):
            plt.plot(days[split_index:], isotope_activities, label=isotopes[idx], linewidth=1)

    plt.plot(days[split_index:], total_activity[split_index:],
             label='Total Activity', linewidth=2, color='black', linestyle='--')

    plt.xlabel('Time (days)')
    plt.ylabel('Activity')
    plt.title(f'{reactor_name} Spent Fuel Activity\nTop 10 Time Integrated Contributors')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.xscale('log')
    plt.yscale('log')
    if late_ylim:
        plt.ylim(late_ylim)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize='small')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_late.png', bbox_inches='tight', dpi=500)
    plt.close()


# Define paths for both reactor types
htgr_base_path = "../../Results/HTGR_FCM/reactor_simulation/"
pwr_base_path = "/home/owen/CLionProjects/R2R4SNF/Results/Ref_PWR/reactor_simulation/"
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

create_activity_plots(htgr_days, htgr_activities, htgr_isotopes,
                      htgr_split_timestep, output_path + 'htgr_activity_plot', 'HTGR FCM',
                      early_ylim=htgr_early_ylim, late_ylim=htgr_late_ylim)
create_activity_plots(pwr_days, pwr_activities, pwr_isotopes,
                      pwr_split_timestep, output_path + 'pwr_activity_plot', 'Reference PWR',
                      early_ylim=pwr_early_ylim, late_ylim=pwr_late_ylim)