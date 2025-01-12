import numpy as np
import matplotlib.pyplot as plt
from plot_helper import *


####################################################################################################
# Common configurations between plots

#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6


####################################################################################################
# Plot 1: B0 vs k0, M1=M2, T=0.0 GeV, mu1=mu2=0.0 GeV, k=0.0 GeV with Mass Shift

# Load data from the file
filename = "B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.0.dat"
data_B0_vs_k0_kvec00 = B0VSK0Data( "../data/" + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
ax.plot(data_B0_vs_k0_kvec00.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec00.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color='red', linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_kvec00.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec00.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color='red', linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
#plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= -2.5; xmax = 2.5; ymin = -5.0; ymax = 5.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=14, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065; auxX = 0.550; auxY = 0.050
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$\mu_{i,j}/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
    r'$|\mathbf{k}|/\Lambda = 0.0$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=14)

auxY=0.95
add_annotation(ax, xmin, xmax, ymin, ymax, auxX, auxY, r'$\mathrm{Shift:}\ M^2 \rightarrow M^2 - i\epsilon$', fontsize=14)

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = filename.replace(".dat", ".png")
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, data_B0_vs_k0_kvec00, plotname
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax


####################################################################################################
# Plot 2: B0 vs k0, M1=M2, T=0.0 GeV, mu1=mu2=0.0 GeV, k=0.0 GeV with with k0 Shift

# Load data from the file
filename = "B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.0_k0Shift.dat"
data_B0_vs_k0_kvec00 = B0VSK0Data( "../data/" + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
ax.plot(data_B0_vs_k0_kvec00.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec00.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color='red', linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_kvec00.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec00.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color='red', linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)

# Configure axes using the helper function
xmin= -2.5
xmax = 2.5
ymin = -5.0
ymax = 5.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=14, spine_width=1.5, tick_width=1.5, tick_length=6)
# Add text annotations
auxH = 0.065
auxX = 0.550
auxY = 0.050
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$\mu_{i,j}/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
    r'$|\mathbf{k}|/\Lambda = 0.0$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=14)

auxY=0.95
add_annotation(ax, xmin, xmax, ymin, ymax, auxX, auxY, r'$\mathrm{Shift:}\ k_0 \rightarrow k_0 + i\epsilon$', fontsize=14)

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = filename.replace(".dat", ".png")
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, data_B0_vs_k0_kvec00, plotname
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax


####################################################################################################
# Plot 3: B0 vs k0, M1=M2, T=0.0 GeV, mu1=mu2=0.0 GeV, k=0.5 GeV with Mass Shift

# Load data from the file
filename = "B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.5.dat"
data_B0_vs_k0_kvec05 = B0VSK0Data( "../data/" + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
ax.plot(data_B0_vs_k0_kvec05.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec05.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color='red', linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_kvec05.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec05.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color='red', linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
# Show legend
#plt.legend(loc='upper left', fontsize=14, frameon=False)

# Configure axes using the helper function
xmin= -2.5
xmax = 2.5
ymin = -5.0
ymax = 5.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=14, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065
auxX = 0.550
auxY = 0.050
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$\mu_{i,j}/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
    r'$|\mathbf{k}|/\Lambda = 0.5$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=14)

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = filename.replace(".dat", ".png")
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, data_B0_vs_k0_kvec05, plotname
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax


####################################################################################################
# Plot 4: B0 vs k0, M1=M2, T=0.0 GeV, mu1=mu2=0.0 GeV, k=1.0 GeV with Mass Shift

# Load data from the file
filename = "B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k1.0.dat"
data_B0_vs_k0_kvec10 = B0VSK0Data( "../data/" + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
ax.plot(data_B0_vs_k0_kvec10.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec10.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color='red', linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_kvec10.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec10.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color='red', linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
# Show legend
#plt.legend(loc='upper left', fontsize=14, frameon=False)

# Configure axes using the helper function
xmin= -2.5
xmax = 2.5
ymin = -5.0
ymax = 5.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=14, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065
auxX = 0.550
auxY = 0.050
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$\mu_{i,j}/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
    r'$|\mathbf{k}|/\Lambda = 1.0$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=14)

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = filename.replace(".dat", ".png")
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, data_B0_vs_k0_kvec10, plotname
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax


####################################################################################################
# Plot 5: B0 vs k0, M1=M2, T=0.0 GeV, mu1=mu2=0.0 GeV, different k with Mass Shift

# Load data from the files
data_B0_vs_k0_kvec00 = B0VSK0Data( "../data/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.0.dat")
data_B0_vs_k0_kvec05 = B0VSK0Data( "../data/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.5.dat")
data_B0_vs_k0_kvec10 = B0VSK0Data( "../data/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k1.0.dat")
data_B0_vs_k0_kvec15 = B0VSK0Data( "../data/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k1.5.dat")

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
ax.plot(data_B0_vs_k0_kvec00.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec00.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color='black', linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_kvec00.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec00.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color='black', linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_kvec05.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec05.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color='red', linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_kvec05.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec05.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color='red', linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_kvec10.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec10.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color='lime', linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_kvec10.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec10.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color='lime', linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_kvec15.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec15.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color='blue', linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_kvec15.get_k0_3D_cutoff_ratio(), 
        data_B0_vs_k0_kvec15.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color='blue', linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
#plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 2.5; ymin = -5.0; ymax = 5.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=14, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065; auxX = 0.050; auxY = 0.050
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$\mu_{i,j}/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)


# Add legend lines and labels
annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, 
                      r'$|\mathbf{k}|/\Lambda = 0.0$', 
                      x_start_factor=0.05, 
                      y_factor=0.94 - auxH*0,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1="black", color2="black", 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, 
                      r'$|\mathbf{k}|/\Lambda = 0.5$', 
                      x_start_factor=0.05, 
                      y_factor=0.94 - auxH*1,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1="red", color2="red", 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, 
                      r'$|\mathbf{k}|/\Lambda = 1.0$', 
                      x_start_factor=0.45, 
                      y_factor=0.94 - auxH*0,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1="lime", color2="lime", 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, 
                      r'$|\mathbf{k}|/\Lambda = 1.5$', 
                      x_start_factor=0.45, 
                      y_factor=0.94 - auxH*1,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1="blue", color2="blue", 
                      width1=2, width2=2, 
                      style1="-", style2="--")

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4_diff_k.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del data_B0_vs_k0_kvec00, data_B0_vs_k0_kvec05, data_B0_vs_k0_kvec10, data_B0_vs_k0_kvec15, plotname
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax
