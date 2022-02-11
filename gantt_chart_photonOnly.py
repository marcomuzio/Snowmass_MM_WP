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

particles = ['gamma']
bands = ['subMeV', 'MeV', 'GeV', 'TeV+']
lgELB = {'subMeV' : -100,
		 'MeV' : 5,
		 'GeV' : 7.5,
		 'TeV+' : 11.5}
lgEUB = {'subMeV' : 5,
		 'MeV' : 7.5,
		 'GeV' : 11.5,
		 'TeV+' : 100}

facecolors = {'subMeV' : 'cornflowerblue',
			  'MeV' : 'red', #'yellow',
			  'GeV' : 'xkcd:gold', #'yellow', #'darkorange',
			  'TeV+' : 'forestgreen'}
edgecolors = {'subMeV' : setAlpha('cornflowerblue',0.8),
			  'MeV' : setAlpha('xkcd:pastel red',0.8),
			  'GeV' : setAlpha('xkcd:gold',0.8),
			  'TeV+' : setAlpha('forestgreen',0.8)}
hatches = {'subMeV' : '///',
		   'MeV' : '\\\\\\',
		   'GeV' : '///',
		   'TeV+' : '\\\\\\'}
bandLbls = {'subMeV' : r'$<$MeV',
			'MeV' : 'MeV',
			'GeV' : 'GeV',
			'TeV+' : r'$\geq$TeV'}


regions = [{},{}]
contributingExps = [{},{}]


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

	# filter out other messengers
	if(not (particle in particles)):
		continue

	lgElow = data['lgElow'].values[i]
	lgEhi = data['lgEhi'].values[i]
	bandRange = [] 
	for band in bands:
		if(max(lgElow, lgELB[band]) <= lgEUB[band] and min(lgEhi, lgEUB[band]) >= lgELB[band]):
			bandRange.append(band)

	isFunded = int(data['isFunded'].values[i])

	Elow = 10**lgElow
	Ehi = 10**lgEhi

	for thisBand in bandRange:

		thisElow = max(Elow, 10**lgELB[thisBand])
		thisEhi = min(Ehi, 10**lgEUB[thisBand])

		r = sg.box(thisElow, start, thisEhi, end)
		if thisBand in regions[isFunded]:
			oldr = regions[isFunded][thisBand]
			regions[isFunded][thisBand] = so.unary_union([oldr,r])
		else:
			regions[isFunded][thisBand] = r

# collect experiment labels for final regions
for i in range(0, len(data)):

	experiment = data['Experiment'].values[i]
	
	start = data['Start Year'].values[i]
	end = data['End Year'].values[i]

	particle = data['Particle'].values[i]

	# filter out other messengers
	if(not (particle in particles)):
		continue

	lgElow = data['lgElow'].values[i]
	lgEhi = data['lgEhi'].values[i]
	bandRange = [] 
	for band in bands:
		if(max(lgElow, lgELB[band]) <= lgEUB[band] and min(lgEhi, lgEUB[band]) >= lgELB[band]):
			bandRange.append(band)

	isFunded = int(data['isFunded'].values[i])

	Elow = 10**lgElow
	Ehi = 10**lgEhi

	for thisBand in bandRange:

		thisElow = max(Elow, 10**lgELB[thisBand])
		thisEhi = min(Ehi, 10**lgEUB[thisBand])

		r = sg.box(thisElow, start, thisEhi, end)
		
		# collect experiments for all but subMeV bands
		if(thisBand == 'subMeV'):
			continue

		# multiple disconnected regions case
		if(regions[isFunded][thisBand].geom_type == 'MultiPolygon'):
			if( not (thisBand in contributingExps[isFunded])):
				contributingExps[isFunded][thisBand] = [[] for j in range(len(regions[isFunded][thisBand].geoms))]

			# check for overlap with each polygon
			for j in range(0, len(regions[isFunded][thisBand].geoms)):
				if(regions[isFunded][thisBand].geoms[j].intersects(r)):
					contributingExps[isFunded][thisBand][j].append(experiment)
		
		# single connected region case
		else:
			if( not (thisBand in contributingExps[isFunded])):
				contributingExps[isFunded][thisBand] = []
			if(regions[isFunded][thisBand].intersects(r)):
				contributingExps[isFunded][thisBand].append(experiment)

# modify unfunded region to remove overlaps with funded region
for band in bands:

	if(not (band in regions[0])):
		continue

	unfunded = regions[0][band]
	funded = regions[1][band]
	
	if(unfunded.intersects(funded)==True):

		nonoverlap = (unfunded.symmetric_difference(funded)).difference(funded)
		regions[0][band] = nonoverlap

print("plotting")
l = [None]*len(bands)
i = 0
for band in bands:

	if(not (band in regions[1])):
		continue

	patch = PolygonPatch(regions[1][band], fc=facecolors[band], ec='dimgrey', alpha=0.25)
	l[i] = plt.gca().add_patch(patch)
	i += 1

for band in bands:
	
	if(not (band in regions[0])):
		continue

	if(regions[0][band].is_empty):
		continue

	patch = PolygonPatch(regions[0][band], fc=(1,1,1,0.3), ec=edgecolors[band], hatch=hatches[band])#, alpha=0.8)
	plt.gca().add_patch(patch)

# make labels and plot them
### no guarantees label placement will always work well ###
for isFunded in range(0, len(contributingExps)):
	for thisBand in contributingExps[isFunded]:

		if(regions[isFunded][thisBand].is_empty):
			continue
		
		# get label coordinates
		if(regions[isFunded][thisBand].geom_type == 'MultiPolygon'):
			expLbl = [None]*len(regions[isFunded][thisBand].geoms)
			x = [None]*len(regions[isFunded][thisBand].geoms)
			y = [None]*len(regions[isFunded][thisBand].geoms)
			for j in range(0, len(regions[isFunded][thisBand].geoms)):
				thisx,thisy = regions[isFunded][thisBand].geoms[j].exterior.coords.xy
				x[j] = np.asarray(thisx).max()/10.**0.1
				y[j] = min(np.asarray(thisy).max(),ymax)
		else:
			expLbl = [None]*1
			x = [None]*1
			y = [None]*1
			thisx,thisy = regions[isFunded][thisBand].exterior.coords.xy
			x[0] = np.asarray(thisx).max()/10.**0.1
			y[0] = min(np.asarray(thisy).max(),ymax)

		# make label
		polyList = contributingExps[isFunded][thisBand]
		if(isinstance(polyList[0],list)): # check for multiple polygons
			for j in range(0, len(polyList)): # loop over polygons if there are multiple
				expLbl[j] = ""
				for k in range(0, len(polyList[j])):
					if(k == 0):
						expLbl[j] += polyList[j][k]
					else:
						expLbl[j] += str(", %s" % polyList[j][k])

		else: # single polygon case
			expLbl[0] = ""
			for j in range(0, len(polyList)):
				if(j == 0):
					expLbl[0] += polyList[j]
				else:
					expLbl[0] += str(", %s" % polyList[j])

		for j in range(0, len(expLbl)):
			if(expLbl[j] is None):
				continue
			if(len(expLbl[j])/13.0 > (y[j]-ymin)/7.0): # try to adjust label if it seems too long (this is definitely not going to work in general)
				buff = expLbl[j].split(', ')
				expLbl[j] = buff[0]
				tmpLen = len(buff[0])
				for k in range(1, len(buff)):
					if((len(buff[k])+tmpLen+2)/13.0 > (y[j]-ymin)/7.0):
						expLbl[j] += ",\n"
						tmpLen = 0
					else:
						expLbl[j] += ", "
						tmpLen += 2
					expLbl[j] += buff[k]
					tmpLen += len(buff[k])
			plt.text(x[j],y[j], expLbl[j], color=facecolors[thisBand], rotation=270, ha='right', va='top')

# add demarcations for EM bands
plt.axvline(x=124e3, linestyle='--', c='grey', lw=1)
plt.axvline(x=124, linestyle='--', c='grey', lw=1)
plt.axvline(x=1.24, linestyle='--', c='grey', lw=1)
plt.axvline(x=1.24e-3, linestyle='--', c='grey', lw=1)
plt.axvline(x=12.4e-6, linestyle='--', c='grey', lw=1)

plt.text(1e-7, 2039.9, 'Radio', rotation=0, fontweight='bold', fontstyle='italic', ha='center', va='top')
plt.text(1.24*10**-1.5, 2039.9, 'IR', rotation=0, fontweight='bold', fontstyle='italic', ha='center', va='top')
plt.text(1.24e1, 2039.9, 'UV', rotation=0, fontweight='bold', fontstyle='italic', ha='center', va='top')
plt.text(124*10**1.5, 2039.9, 'X-ray', rotation=0, fontweight='bold', fontstyle='italic', ha='center', va='top')

#plt.text(1.24e-3, 2039.5, 'Infrared', rotation=270, fontweight='bold', fontstyle='italic', ha='left', va='top')
#plt.text(1.24, 2039.5, 'UV', rotation=270, fontweight='bold', fontstyle='italic', ha='left', va='top')
#plt.text(124, 2039.5, 'X-ray', rotation=270, fontweight='bold', fontstyle='italic', ha='left', va='top')

plt.xlim(1e-9, 1e22)
plt.ylim(ymin, ymax)

plt.xscale('log')
plt.gca().xaxis.set_ticks(10.**np.arange(-9,22,3))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%g'))
plt.gca().yaxis.set_ticks(range(ymin, ymax+1, 5))

plt.xlabel('Energy (eV)')
plt.ylabel('Year')


new_tick_locations = 10.**np.arange(6,37,3)
ax1 = plt.gca()
ax2 = ax1.twiny()
ax2.set_xlim(np.array(ax1.get_xlim())/hPlanck)
ax2.set_xscale('log')
ax2.set_xticks(new_tick_locations)
ax2.set_xlabel(r"Frequency (Hz)")


# fake line for legend
lbls = [None]*len(bands)
i = 0
for band in bands:
	lbls[i] = bandLbls[band]
	i += 1

plt.legend(handles=l,labels=lbls, loc=8, bbox_to_anchor=(0.5, 1.1, 0., 0.), numpoints=1, ncol=4, frameon=False, framealpha=0.)

plt.tight_layout()

plt.draw()
plt.show()
