import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set up LaTeX fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

# Style settings
background_color = '#FFFFFF'
grid_color = '#333333'
text_color = '#333333'

# Font sizes
TITLE_SIZE = 16
LABEL_SIZE = 14
TICK_SIZE = 12
LEGEND_SIZE = 12

# Read the CSV data
df = pd.read_csv('../../Results/I129_release_transport/PFLOTRAN_output.csv')

# Create figure
fig, ax = plt.subplots(figsize=(12, 7))

# Plot each reactor type
reactors = ['Ref_PWR', 'HTGR_FCM']
colors = ['#000000', '#000000', '#2ca02c', '#d62728', '#9467bd', '#8c564b']  # More professional color palette
styles = ['-', '--', ':', '-.', '-', '--']

for reactor, color, style in zip(reactors, colors, styles):
    plt.semilogy(df['Time [y]'], df[reactor], label=reactor, 
                 color=color, linestyle=style, linewidth=1.5)

# Style the plot
ax.set_facecolor(background_color)
fig.patch.set_facecolor(background_color)
ax.grid(True, which="both", ls="-", alpha=0.2, color=grid_color)
ax.set_axisbelow(True)

# Style spines
for spine in ax.spines.values():
    spine.set_color('#CCCCCC')

# Style labels and ticks
ax.tick_params(colors=text_color, labelsize=TICK_SIZE)
ax.set_xlabel('Time (years)', fontsize=LABEL_SIZE, color=text_color)
ax.set_ylabel('Concentration', fontsize=LABEL_SIZE, color=text_color)

# Create and style legend
leg = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                borderaxespad=0., fontsize=LEGEND_SIZE, framealpha=1)
leg.get_frame().set_facecolor(background_color)

output_path = "../output/"

# Adjust layout and save
plt.tight_layout()
plt.savefig(output_path + 'pflotran_comparison.png', bbox_inches='tight', dpi=500,
            facecolor=background_color, edgecolor='none')
plt.close()