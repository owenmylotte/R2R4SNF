import matplotlib.pyplot as plt
import numpy as np

# %% Plot radioactivity from SNF based on values in ../../Results/HTGR_FCM/reactor_simulation/

activity_directory_path = "../../Results/HTGR_FCM/reactor_simulation/"
output_path = "../output/"

# The below functions only need to exist because the data output format is ... non-standard.
def read_horizontal_data(filename):
    """
    Reads horizontally oriented data where each row represents a different measurement
    and columns are different parameters.

    Returns:
    - header: list of column names
    - days: numpy array of days
    - data: numpy array of data values
    """
    with open(filename, 'r') as f:
        # Read header line
        header = f.readline().strip().split(',')

        # Read the rest into a list of lists
        data_rows = []
        for line in f:
            # Convert all values to float
            try:
                row = [float(val) for val in line.strip().split(',')]
                data_rows.append(row)
            except ValueError:
                continue

        # Convert to numpy array
        data = np.array(data_rows)

        # Extract days column (first column)
        days = data[:, 0]

    return header, days, data

def read_vertical_data(filename):
    """
    Reads vertically oriented data where first column contains labels
    and subsequent columns contain time series data.

    Returns:
    - time_points: list of column headers (time points)
    - isotopes: list of isotope names
    - activities: numpy array of activity values
    """
    with open(filename, 'r') as f:
        # Read header line for time points
        time_points = f.readline().strip().split(',')

        # Initialize lists for isotopes and data
        isotopes = []
        data_rows = []

        # Read each line
        for line in f:
            values = line.strip().split(',')
            isotope = values[0]

            # Convert remaining values to float where possible
            try:
                row_data = [float(val) if val.strip() else np.nan for val in values[1:]]
                isotopes.append(isotope)
                data_rows.append(row_data)
            except ValueError:
                continue

        # Convert to numpy array
        activities = np.array(data_rows)

    return time_points, isotopes, activities


# Read the data
days_header, days, days_data = read_horizontal_data(activity_directory_path + 'summary/out_HTGR_FCM_time_dep.csv')
time_points, isotopes, activities = read_vertical_data(activity_directory_path + 'SNF_by_nuclide/radioactivity/HTGR_FCM_activity.csv')

# Create figure and axis
plt.figure(figsize=(12, 8))

# Plot each isotope's activity over time
for i, isotope in enumerate(isotopes):
    # Get activities for this isotope
    isotope_activities = activities[i, :]

    # Only plot if we have valid numerical data
    valid_indices = ~np.isnan(isotope_activities)
    if np.any(valid_indices):
        valid_days = days[:np.sum(valid_indices)]
        valid_activities = isotope_activities[valid_indices]

        if len(valid_days) > 0 and len(valid_activities) > 0:
            plt.plot(valid_days, valid_activities, label=isotope, linewidth=1)


# TODO: Maybe we should limit the plotting or legend to only the most prominent nuclides.
#  The question is open ended enough that it makes sense to include only the top contributors.
#  It may also make sense to split the plot into one of each for before and after shutdown.

# Customize the plot
plt.xlabel('Time (days)')
plt.ylabel('Activity')
plt.title('Nuclear Activity Over Time for Different Isotopes')
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.xscale('log')  # Use log scale for better visualization of wide range of values

# Adjust legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize='small')

# Adjust layout to prevent legend cutoff
plt.tight_layout()

# Make the output directory if it does not exist
import os
os.makedirs(output_path, exist_ok=True)

# Save the plot
plt.savefig(output_path + 'activity_plot.png', bbox_inches='tight', dpi=500)
plt.close()

# %% Plot fissile concentrations by nuclide.

# TODO: This will be essentially the same as the previous plot, but with different data for concentrations.
#  The same parsing functions can be used.
