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
# B0 vs k0, M1=M2, T=0.0 GeV, k=0.0 GeV, different mu1=mu2 with Mass Shift
print("Building plot: B0 vs k0, M1=M2, T=0.0 GeV, k=0.0 GeV, different mu1=mu2 with Mass Shift")

# Load data from the files
data_B0_vs_k0_mu00 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.0.dat")
data_B0_vs_k0_mu06 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.6Cpj0.6L1.0Mi0.4Mj0.4k0.0.dat")
data_B0_vs_k0_mu08 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.8Cpj0.8L1.0Mi0.4Mj0.4k0.0.dat")
data_B0_vs_k0_mu10 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi1.0Cpj1.0L1.0Mi0.4Mj0.4k0.0.dat")

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
colorT00 = 'black'
colorT03 = 'red'
colorT05 = 'lime'
colorT10 = 'blue'

ax.plot(data_B0_vs_k0_mu00.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu00.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT00, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu00.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu00.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT00, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu06.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu06.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT03, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu06.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu06.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT03, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu08.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu08.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT05, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu08.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu08.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT05, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu10.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu10.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT10, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu10.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu10.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT10, linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
#plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 2.5; ymin = -4.0; ymax = 6.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=6, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065; auxX = 0.030; auxY = 0.74
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
    r'$|\mathbf{k}|/\Lambda = 0.0$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)


# Add legend lines and labels
dist = 0.06

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.0$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*0,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT00, color2=colorT00, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.6$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*1,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT03, color2=colorT03, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.8$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*2,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT05, color2=colorT05, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 1.0$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*3,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT10, color2=colorT10, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "B0_vs_k0_T0.0L1.0Mi0.4Mj0.4k0.0_diff_mu.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del data_B0_vs_k0_mu00, data_B0_vs_k0_mu06, data_B0_vs_k0_mu08, data_B0_vs_k0_mu10, plotname, dist
del colorT00, colorT03, colorT05, colorT10
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax


####################################################################################################
# B0 vs k0, M1=M2, T=0.0 GeV, k=0.5 GeV, different mu1=mu2 with Mass Shift
print("Building plot: B0 vs k0, M1=M2, T=0.0 GeV, k=0.5 GeV, different mu1=mu2 with Mass Shift")

# Load data from the files
data_B0_vs_k0_mu00 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.5.dat")
data_B0_vs_k0_mu06 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.6Cpj0.6L1.0Mi0.4Mj0.4k0.5.dat")
data_B0_vs_k0_mu08 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.8Cpj0.8L1.0Mi0.4Mj0.4k0.5.dat")
data_B0_vs_k0_mu10 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi1.0Cpj1.0L1.0Mi0.4Mj0.4k0.5.dat")

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
colorT00 = 'black'
colorT03 = 'red'
colorT05 = 'lime'
colorT10 = 'blue'

ax.plot(data_B0_vs_k0_mu00.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu00.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT00, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu00.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu00.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT00, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu06.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu06.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT03, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu06.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu06.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT03, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu08.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu08.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT05, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu08.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu08.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT05, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu10.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu10.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT10, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu10.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu10.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT10, linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
#plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 2.5; ymin = -4.0; ymax = 6.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=6, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065; auxX = 0.030; auxY = 0.74
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
    r'$|\mathbf{k}|/\Lambda = 0.5$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)


# Add legend lines and labels
dist = 0.06

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.0$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*0,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT00, color2=colorT00, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.6$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*1,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT03, color2=colorT03, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.8$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*2,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT05, color2=colorT05, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 1.0$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*3,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT10, color2=colorT10, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "B0_vs_k0_T0.0L1.0Mi0.4Mj0.4k0.5_diff_mu.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del data_B0_vs_k0_mu00, data_B0_vs_k0_mu06, data_B0_vs_k0_mu08, data_B0_vs_k0_mu10, plotname, dist
del colorT00, colorT03, colorT05, colorT10
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax


####################################################################################################
# B0 vs k0, M1=M2, T=0.0 GeV, k=1.0 GeV, different mu1=mu2 with Mass Shift
print("Building plot: B0 vs k0, M1=M2, T=0.0 GeV, k=1.0 GeV, different mu1=mu2 with Mass Shift")

# Load data from the files
data_B0_vs_k0_mu00 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k1.0.dat")
data_B0_vs_k0_mu06 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.6Cpj0.6L1.0Mi0.4Mj0.4k1.0.dat")
data_B0_vs_k0_mu08 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.8Cpj0.8L1.0Mi0.4Mj0.4k1.0.dat")
data_B0_vs_k0_mu10 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi1.0Cpj1.0L1.0Mi0.4Mj0.4k1.0.dat")

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
colorT00 = 'black'
colorT03 = 'red'
colorT05 = 'lime'
colorT10 = 'blue'

ax.plot(data_B0_vs_k0_mu00.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu00.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT00, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu00.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu00.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT00, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu06.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu06.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT03, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu06.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu06.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT03, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu08.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu08.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT05, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu08.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu08.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT05, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu10.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu10.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT10, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu10.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu10.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT10, linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
#plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 2.5; ymin = -4.0; ymax = 6.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=6, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065; auxX = 0.030; auxY = 0.74
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
    r'$|\mathbf{k}|/\Lambda = 1.0$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)


# Add legend lines and labels
dist = 0.06

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.0$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*0,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT00, color2=colorT00, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.6$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*1,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT03, color2=colorT03, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.8$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*2,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT05, color2=colorT05, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 1.0$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*3,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT10, color2=colorT10, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "B0_vs_k0_T0.0L1.0Mi0.4Mj0.4k1.0_diff_mu.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del data_B0_vs_k0_mu00, data_B0_vs_k0_mu06, data_B0_vs_k0_mu08, data_B0_vs_k0_mu10, plotname, dist
del colorT00, colorT03, colorT05, colorT10
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax


####################################################################################################
# B0 vs k0, M1=M2, T=0.0 GeV, k=1.5 GeV, different mu1=mu2 with Mass Shift
print("Building plot: B0 vs k0, M1=M2, T=0.0 GeV, k=1.5 GeV, different mu1=mu2 with Mass Shift")

# Load data from the files
data_B0_vs_k0_mu00 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k1.5.dat")
data_B0_vs_k0_mu06 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.6Cpj0.6L1.0Mi0.4Mj0.4k1.5.dat")
data_B0_vs_k0_mu08 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi0.8Cpj0.8L1.0Mi0.4Mj0.4k1.5.dat")
data_B0_vs_k0_mu10 = B03DCutoffVsMomentumData( "../data/B0_vs_k0_T0.0Cpi1.0Cpj1.0L1.0Mi0.4Mj0.4k1.5.dat")

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
colorT00 = 'black'
colorT03 = 'red'
colorT05 = 'lime'
colorT10 = 'blue'

ax.plot(data_B0_vs_k0_mu00.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu00.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT00, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu00.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu00.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT00, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu06.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu06.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT03, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu06.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu06.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT03, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu08.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu08.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT05, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu08.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu08.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT05, linewidth=2, linestyle='--')

ax.plot(data_B0_vs_k0_mu10.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu10.get_Re_B0(), 
        label=r'$\mathrm{Re}[B_0]$', color=colorT10, linewidth=2, linestyle='-')
ax.plot(data_B0_vs_k0_mu10.get_momentum_to_Lambda_ratio(), 
        data_B0_vs_k0_mu10.get_Im_B0(), 
        label=r'$\mathrm{Im}[B_0]$', color=colorT10, linewidth=2, linestyle='--')

# Axes labels
ax.set_xlabel(r'$k_0/\Lambda$', fontsize=20)
ax.set_ylabel(r'$B_0$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
#plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 2.5; ymin = -4.0; ymax = 6.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=6, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

# Add text annotations
auxH = 0.065; auxX = 0.030; auxY = 0.74
texts = [
    r'$\Lambda = 1.0\ \mathrm{GeV}$',
    r'$T/\Lambda = 0.0$',
    r'$M_{i,j}/\Lambda = 0.4$',
    r'$|\mathbf{k}|/\Lambda = 1.5$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)


# Add legend lines and labels
dist = 0.06

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.0$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*0,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT00, color2=colorT00, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.6$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*1,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT03, color2=colorT03, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 0.8$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*2,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT05, color2=colorT05, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, dist,
                      r'$\mu_{i,j}/\Lambda = 1.0$', 
                      x_start_factor=0.03, 
                      y_factor=0.24 - auxH*3,
                      length_line_fraction_of_box=0.4, 
                      fontsize=16, 
                      color1=colorT10, color2=colorT10, 
                      width1=2, width2=2, 
                      style1="-", style2="--")

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "B0_vs_k0_T0.0L1.0Mi0.4Mj0.4k1.5_diff_mu.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del data_B0_vs_k0_mu00, data_B0_vs_k0_mu06, data_B0_vs_k0_mu08, data_B0_vs_k0_mu10, plotname, dist
del colorT00, colorT03, colorT05, colorT10
del fig, ax, auxH, auxX, auxY, xmin, xmax, ymin, ymax
