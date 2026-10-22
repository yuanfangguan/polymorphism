import matplotlib.pyplot as plt
import numpy as np

# Example data (replace with your actual data)
data=[]
FILE=open("../base_excludeN_5fold/auc.perposition",'r')
for line in FILE:
    line=line.strip()
    data.append(float(line))
FILE.close()
# Plot the first 20 and last 20 positions with a broken axis in the middle
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3), gridspec_kw={'width_ratios': [1, 1]})

# First 20 positions
ax1.plot(range(50), data[:50], color='blue')
ax1.set_xlim(0, 49)
ax1.set_title('First 50 Positions')

# Last 20 positions
ax2.plot(range(1998, 2048), data[-50:], color='blue')
ax2.set_xlim(1998, 2047)
ax2.set_title('Last 50 Positions')

# Hide the spines between ax1 and ax2
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
ax2.yaxis.tick_right()

# Add diagonal lines to indicate the break
d = 0.025  # how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d, 1+d), (-d, +d), **kwargs)        # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal

# Labels and titles
fig.text(0.5, 0.04, 'Position', ha='center')
fig.text(0.04, 0.5, 'AUC', va='center', rotation='vertical')
#fig.suptitle('Line Plot with Expanded Tails and Compressed Middle')

# Save the figure as a PNG file
plt.savefig('auc_position.png', dpi=300)

# Show the plot
plt.show()

