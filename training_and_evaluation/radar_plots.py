import numpy as np
import matplotlib.pyplot as plt

# =========================
# DATA 
# =========================
categories = ["Mixed", "Drug-Blind", "Cancer-Blind", "Disjoint"]

gdrp_old = [.75, -.11, .60, -.26] # add results for plot, in this order : ["Mixed", "Drug-Blind", "Cancer-Blind", "Disjoint"]
gdrp_new = [.66, .12, .57, .08]

# =========================
# CONFIGURATION
# =========================
r_min = -0.4
r_max = 1.0
r_step = 0.2

title = "GraphDRP 10-Fold Mean R²"

# styling
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.labelsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 10,
    "font.family": "DejaVu Sans",
})

# =========================
# RADAR SETUP
# =========================
N = len(categories)
angles = np.linspace(0, 2 * np.pi, N, endpoint=False)

# Close the loop
angles = np.concatenate((angles, [angles[0]]))
gdrp_old = np.concatenate((gdrp_old, [gdrp_old[0]]))
gdrp_new = np.concatenate((gdrp_new, [gdrp_new[0]]))

fig, ax = plt.subplots(figsize=(8, 7), subplot_kw=dict(polar=True))

# Rotate so first axis is at top
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# =========================
# AXIS & GRID
# =========================
ax.set_ylim(r_min, r_max)

radial_ticks = np.arange(r_min, r_max + r_step, r_step)
ax.set_yticks(radial_ticks)


ax.set_yticklabels(
    ["" if tick == r_min else f"{tick:.1f}" for tick in radial_ticks]
)

ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories)

# Move specific labels outward
for label, angle in zip(ax.get_xticklabels(), angles[:-1]):
    if label.get_text() in ["Drug-Blind", "Disjoint"]:
        label.set_y(label.get_position()[1] - 0.13)

# grid
# ax.grid(color="black", linestyle="--", linewidth=0.6, alpha=0.8)
ax.yaxis.grid(True, color="black", linestyle="--", linewidth=0.6, alpha=0.8)
ax.xaxis.grid(True, color="gray", linestyle="--", linewidth=0.6, alpha=0.6)

ax.spines["polar"].set_visible(False)

# =========================
# PLOTTING
# =========================
color_original = "#d62728"
color_expanded = "#1f77b4"

ax.plot(angles, gdrp_old, linewidth=2.5, label="Original", color=color_original)
ax.fill(angles, gdrp_old, alpha=0.15, color=color_original)

ax.plot(angles, gdrp_new, linewidth=2.5, label="Expanded", color=color_expanded)
ax.fill(angles, gdrp_new, alpha=0.15, color=color_expanded)

# =========================
# TITLE & LEGEND
# =========================
ax.set_title(title, pad=45)
ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1.10), frameon=False)

plt.tight_layout()
plt.savefig(
    "gdrp_radarplot.png",
    dpi=300
)
plt.show()
