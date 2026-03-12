import numpy as np
import matplotlib.pyplot as plt
from common_utils.plot_helper import *
from common_utils.b03d_cutoff_vs_momentum_data import *


####################################################################################################
# Common configurations between plots

#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
data_folder = "two_fermion_line_integral_3d_cutoff/data/"
plots_folder = "two_fermion_line_integral_3d_cutoff/plots/"

####################################################################################################
# B0 vs k, M1=M2, T=0.0 GeV, mu1=mu2=0.0 GeV, different k0 with Mass Shift
print("Building plot: B0 vs k, M1=M2, T=0.0 GeV, mu1=mu2=0.0 GeV, different k0 with Mass Shift")

# Load data from the file
data_B0_vs_k_k000 = B03DCutoffVsMomentumData(data_folder + "B0_vs_k_T0p0Cpi0p0Cpj0p0L1p0Mi0p4Mj0p4k00p0.dat")
data_B0_vs_k_k005 = B03DCutoffVsMomentumData(data_folder + "B0_vs_k_T0p0Cpi0p0Cpj0p0L1p0Mi0p4Mj0p4k00p5.dat")
data_B0_vs_k_k010 = B03DCutoffVsMomentumData(data_folder + "B0_vs_k_T0p0Cpi0p0Cpj0p0L1p0Mi0p4Mj0p4k01p0.dat")
data_B0_vs_k_k015 = B03DCutoffVsMomentumData(data_folder + "B0_vs_k_T0p0Cpi0p0Cpj0p0L1p0Mi0p4Mj0p4k01p5.dat")
data_B0_vs_k_k020 = B03DCutoffVsMomentumData(data_folder + "B0_vs_k_T0p0Cpi0p0Cpj0p0L1p0Mi0p4Mj0p4k02p0.dat")
data_B0_vs_k_k025 = B03DCutoffVsMomentumData(data_folder + "B0_vs_k_T0p0Cpi0p0Cpj0p0L1p0Mi0p4Mj0p4k02p5.dat")

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
colork000 = 'black'
colork005 = 'red'
colork010 = 'lime'
colork015 = 'blue'
colork020 = (0, 0.9, 0.9)  # Light cyan/turquoise
colork025 = (1, 0.1, 0.9)   # Magenta

ax.plot(data_B0_vs_k_k000.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k000.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colork000, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k_k000.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k000.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colork000, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k_k005.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k005.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colork005, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k_k005.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k005.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colork005, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k_k010.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k010.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colork010, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k_k010.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k010.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colork010, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k_k015.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k015.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colork015, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k_k015.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k015.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colork015, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k_k020.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k020.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colork020, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k_k020.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k020.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colork020, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k_k025.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k025.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colork025, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k_k025.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k_k025.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colork025, linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$|\mathbf{k}|/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)

# Configure axes using the helper function
xmin= 0.0; xmax = 2.5; ymin = -3.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=6, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065; auxX = 0.66; auxY = 0.75
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$\mu_{i,j}/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

# Add legend lines and labels
dist = 0.06

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$k_0/\Lambda = 0.0$', 
                      x_start_factor=0.58, 
                      y_factor=0.35 - auxH*0,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colork000, color2=colork000, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$k_0/\Lambda = 0.5$', 
                      x_start_factor=0.58, 
                      y_factor=0.35 - auxH*1,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colork005, color2=colork005, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$k_0/\Lambda = 1.0$', 
                      x_start_factor=0.58, 
                      y_factor=0.35 - auxH*2,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colork010, color2=colork010, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$k_0/\Lambda = 1.5$', 
                      x_start_factor=0.58, 
                      y_factor=0.35 - auxH*3,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colork015, color2=colork015, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$k_0/\Lambda = 2.0$', 
                      x_start_factor=0.58, 
                      y_factor=0.35 - auxH*4,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colork020, color2=colork020, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$k_0/\Lambda = 2.5$', 
                      x_start_factor=0.58, 
                      y_factor=0.35 - auxH*5,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colork025, color2=colork025, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "B0_vs_k_T0p0Cpi0p0Cpj0p0L1p0Mi0p4Mj0p4_diff_k0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del data_B0_vs_k_k000, plotname, dist
del colork000, colork005, colork010, colork015, colork020, colork025
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax
