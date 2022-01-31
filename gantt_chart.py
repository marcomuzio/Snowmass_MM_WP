import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from descartes import PolygonPatch
import shapely.geometry as sg
import shapely.ops as so

data = pd.read_csv("experiments.txt", index_col=False, comment='#')

particles = ['gamma', 'cr', 'nu', 'gw']

colors = {'gamma' : 'cornflowerblue',
		  'cr' : 'yellow',
		  'nu' : 'darkorange',
		  'gw' : 'forestgreen'}

particleLbls = ['Photons', 'CRs', 'Neutrinos', 'GWs']

regions = {}

#cmaps = ['Greens', 'Blues', 'Greys', 'Purples']

hPlanck = 4.136e-15 # eV/Hz

ymin = 2020
ymax = 2040

matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = "stix"

for i in range(0, len(data)):

	experiment = data['Experiment'].values[i]
	
	start = data['Start Year'].values[i]
	end = data['End Year'].values[i]

	particle = data['Particle'].values[i]

	lgElow = data['lgElow'].values[i]
	lgEhi = data['lgEhi'].values[i]

	Elow = 10**lgElow
	Ehi = 10**lgEhi

	yposmin = max(ymin,start)
	yposmax = min(ymax,end)
	ypos = yposmin + 0.75*(yposmax-yposmin)

	#plt.fill_between([Elow, Ehi], start, end, facecolor=colors[particle], edgecolor='dimgrey', label=experiment, alpha=0.5)
	
	print("box", experiment, Elow, start, Ehi, end)
	r = sg.box(Elow, start, Ehi, end)
	if particle in regions:
		oldr = regions[particle]
		regions[particle] = so.unary_union([oldr,r])
	else:
		regions[particle] = r
	#plt.text(Elow*10**0.1, ypos, experiment, color='white', rotation=90)

print("plotting")
for particle in particles:

	patch = PolygonPatch(regions[particle], fc=colors[particle], ec='dimgrey', alpha=0.5)
	plt.gca().add_patch(patch)
	#xs, ys = regions[particle].exterior.coords.xy
	#plt.gca().fill(xs, ys, facecolor=colors[particle], edgecolor='dimgrey', alpha=0.5)

plt.xlim(1e-21, 1e22)
plt.ylim(ymin, ymax)

plt.xscale('log')
plt.gca().xaxis.set_ticks(10.**np.arange(-21,22,3))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%g'))
plt.gca().yaxis.set_ticks(range(ymin, ymax+1, 5))

plt.xlabel('Energy (eV)')
plt.ylabel('Year')


new_tick_locations = 10.**np.arange(-6,37,3)
ax1 = plt.gca()
ax2 = ax1.twiny()
ax2.set_xlim(np.array(ax1.get_xlim())/hPlanck)
ax2.set_xscale('log')
ax2.set_xticks(new_tick_locations)
#ax2.set_xticklabels(new_tick_locations)
ax2.set_xlabel(r"Frequency (Hz)")


# fake line for legend
l = [None]*len(particles)
i = 0
for particle in particles:
	l[i] = plt.fill_between([1e0,1e0], 0, 0, facecolor=colors[particle], edgecolor='dimgrey', alpha=0.5)
	i += 1

plt.legend(handles=l,labels=particleLbls, loc=8, bbox_to_anchor=(0.5, 1.1, 0., 0.), numpoints=1, ncol=4, frameon=False, framealpha=0.) # prop={'size':16},

plt.tight_layout()

plt.draw()
plt.show()
