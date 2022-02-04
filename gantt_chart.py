import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from descartes import PolygonPatch
import shapely.geometry as sg
import shapely.ops as so
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 1.0

def setAlpha(c, alpha):
	c = mpl.colors.to_rgba_array(c)[0]
	c[-1] = alpha
	return c

data = pd.read_csv("experiments.txt", index_col=False, comment='#')

particles = ['gamma', 'nu', 'cr', 'gw']

facecolors = {'gamma' : 'cornflowerblue',
			  'cr' : 'red', #'yellow',
			  'nu' : 'xkcd:gold', #'yellow', #'darkorange',
			  'gw' : 'forestgreen'}
edgecolors = {'gamma' : setAlpha('cornflowerblue',0.8),
			  'cr' : setAlpha('xkcd:pastel red',0.8),
			  'nu' : setAlpha('xkcd:gold',0.8),
			  'gw' : setAlpha('forestgreen',0.8)}
hatches = {'gamma' : '///',
		   'cr' : '\\\\\\',
		   'nu' : '///',
		   'gw' : '\\\\\\'}
particleLbls = {'gamma' : 'Photons',
				'nu' : 'Neutrinos',
				'cr' : 'CRs',
				'gw' : 'GWs'}


regions = [{},{}]


hPlanck = 4.136e-15 # eV/Hz

ymin = 2020
ymax = 2040

matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = "stix"

# get regions
for i in range(0, len(data)):

	experiment = data['Experiment'].values[i]
	
	start = data['Start Year'].values[i]
	end = data['End Year'].values[i]

	particle = data['Particle'].values[i]

	lgElow = data['lgElow'].values[i]
	lgEhi = data['lgEhi'].values[i]

	isFunded = int(data['isFunded'].values[i])

	Elow = 10**lgElow
	Ehi = 10**lgEhi

	yposmin = max(ymin,start)
	yposmax = min(ymax,end)
	ypos = yposmin + 0.75*(yposmax-yposmin)

	#plt.fill_between([Elow, Ehi], start, end, facecolor=colors[particle], edgecolor='dimgrey', label=experiment, alpha=0.5)
	
	#print("box", experiment, Elow, start, Ehi, end)
	r = sg.box(Elow, start, Ehi, end)
	if particle in regions[isFunded]:
		oldr = regions[isFunded][particle]
		regions[isFunded][particle] = so.unary_union([oldr,r])
	else:
		regions[isFunded][particle] = r
	#plt.text(Elow*10**0.1, ypos, experiment, color='white', rotation=90)

# modify unfunded region to remove overlaps with funded region
for particle in particles:

	unfunded = regions[0][particle]
	funded = regions[1][particle]
	
	if(unfunded.intersects(funded)==True):

		nonoverlap = (unfunded.symmetric_difference(funded)).difference(funded)
		regions[0][particle] = nonoverlap

print("plotting")
l = [None]*len(particles)
i = 0
for particle in particles:

	#patch = PolygonPatch(regions[0][particle], fc='white', ec=colors[particle], hatch=hatches[particle], alpha=0.25)
	#patch = PolygonPatch(regions[0][particle], fc=colors[particle], ec='dimgrey', hatch='/', alpha=0.25)
	#plt.gca().add_patch(patch)
	patch = PolygonPatch(regions[1][particle], fc=facecolors[particle], ec='dimgrey', alpha=0.25)
	#patch = PolygonPatch(regions[1][particle], fc='white', ec=colors[particle], hatch=hatches[particle], alpha=0.4)
	l[i] = plt.gca().add_patch(patch)
	i += 1

for particle in particles:

	if(regions[0][particle].is_empty):
		continue

	patch = PolygonPatch(regions[0][particle], fc=(1,1,1,0.3), ec=edgecolors[particle], hatch=hatches[particle])#, alpha=0.8)
	#patch = PolygonPatch(regions[0][particle], fc=colors[particle], ec='dimgrey', hatch='/', alpha=0.25)
	plt.gca().add_patch(patch)
	#patch = PolygonPatch(regions[1][particle], fc=colors[particle], ec='dimgrey', alpha=0.4)
	#plt.gca().add_patch(patch)

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
ax2.set_xlabel(r"Frequency (Hz)")


# fake line for legend
#l = [None]*len(particles)
lbls = [None]*len(particles)
i = 0
for particle in particles:
	#l[i] = plt.fill_between([1e0,1e0], 0, 0, facecolor=facecolors[particle], edgecolor='dimgrey', alpha=0.5)
	lbls[i] = particleLbls[particle]
	i += 1

plt.legend(handles=l,labels=lbls, loc=8, bbox_to_anchor=(0.5, 1.1, 0., 0.), numpoints=1, ncol=4, frameon=False, framealpha=0.)

plt.tight_layout()

plt.draw()
plt.show()
