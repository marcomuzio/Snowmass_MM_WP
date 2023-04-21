import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
from descartes import PolygonPatch
import shapely.geometry as sg
import shapely.ops as so
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 1.0
import sys

def setAlpha(c, alpha):
	c = mpl.colors.to_rgba_array(c)[0]
	c[-1] = alpha
	return c

def wrapwidth(x):
	return lambda : x

data = pd.read_csv("experiments.txt", index_col=False, comment='#')

particles = ['gw', 'gamma', 'nu', 'cr']

textcolors = {'gamma' : 'royalblue',
			  'cr' : 'darkred', #'yellow',
			  'nu' : 'xkcd:dark gold', #'yellow', #'darkorange',
			  'gw' : 'green'}
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

bands = {'gamma' : ['subMeV', 'MeV', 'GeV', 'TeV+'],
		 'cr' : ['subPeV', 'PeV', 'PeV+'],
		 'nu' : ['subPeV', 'PeV', 'PeV+'],
		 'gw' : ['all']}
lgELB = {'subMeV' : -100, 
		 'MeV' : 5,
		 'GeV' : 7.5,
		 'TeV+' : 11.5,
		 'subPeV': -100,
		 'PeV': 14,
		 'PeV+':17,
		 'all' : -100}
lgEUB = {'subMeV' : 5,
		 'MeV' : 7.5,
		 'GeV' : 11.5,
		 'TeV+' : 100,
		 'subPeV': 14,
		 'PeV': 17,
		 'PeV+': 100,
		 'all' : 100}

regions = [{},{}]
contributingExps = [{},{},{}]
axOut = {'gamma' : 0,
		 'cr' : 1,
		 'nu' : 1,
		 'gw' : 1}
axIn = {'gamma' : 1,
		'cr' : 2,
		'nu' : 1,
		'gw' : 0}
	
widths = np.zeros([2,3]) # should match number of panels in gridspec

hPlanck = 4.136e-15 # eV/Hz

ymin = 2020
ymax = 2040

matplotlib.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = "Helvetica"

# get regions
for i in range(0, len(data)):

	experiment = data['Experiment'].values[i]
	
	start = data['Start Year'].values[i]
	if(start < ymin):
		start = ymin
	end = data['End Year'].values[i]
	if(end > ymax):
		end = ymax

	particle = data['Particle'].values[i]

	lgElow = data['lgElow'].values[i]
	lgEhi = data['lgEhi'].values[i]
	bandRange = [] 
	for band in bands[particle]:
		if(max(lgElow, lgELB[band]) <= lgEUB[band] and min(lgEhi, lgEUB[band]) >= lgELB[band]):
			bandRange.append(band)

	isFunded = int(data['isFunded'].values[i])

	Elow = 10**lgElow
	Ehi = 10**lgEhi

	#print("box", experiment, Elow, start, Ehi, end)
	for thisBand in bandRange:

		thisElow = max(Elow, 10**lgELB[thisBand])
		thisEhi = min(Ehi, 10**lgEUB[thisBand])

		r = sg.box(thisElow, start, thisEhi, end)
		if( not (particle in regions[isFunded])):
			regions[isFunded][particle] = {}
		if thisBand in regions[isFunded][particle]:
			oldr = regions[isFunded][particle][thisBand]
			regions[isFunded][particle][thisBand] = so.unary_union([oldr,r])
		else:
			regions[isFunded][particle][thisBand] = r
	
	#r = sg.box(Elow, start, Ehi, end)
	#if particle in regions[isFunded]:
	#	oldr = regions[isFunded][particle]
	#	regions[isFunded][particle] = so.unary_union([oldr,r])
	#else:
	#	regions[isFunded][particle] = r

# modify unfunded region to remove overlaps with funded region
for particle in particles:

	if(not (particle in regions[0])):
		continue

	for band in regions[0][particle]:
		unfunded = regions[0][particle][band]
		funded = regions[1][particle][band]
		
		if(unfunded.intersects(funded)==True):

			nonoverlap = (unfunded.symmetric_difference(funded)).difference(funded)
			regions[0][particle][band] = nonoverlap

# collect experiment labels for final regions
for i in range(0, len(data)):

	experiment = data['Experiment'].values[i]
	
	start = data['Start Year'].values[i]
	if(start < ymin):
		start = ymin
	end = data['End Year'].values[i]
	if(end > ymax):
		end = ymax

	thisParticle = data['Particle'].values[i]

	lgElow = data['lgElow'].values[i]
	lgEhi = data['lgEhi'].values[i]
	bandRange = [] 
	for band in bands[thisParticle]:
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
		#if(thisBand == 'subMeV'):
		#	continue

		# multiple disconnected regions case
		if(regions[isFunded][thisParticle][thisBand].geom_type == 'MultiPolygon'):
			if( not (thisParticle in contributingExps[isFunded])):
				contributingExps[isFunded][thisParticle] = {}
			if( not (thisBand in contributingExps[isFunded][thisParticle])):
				contributingExps[isFunded][thisParticle][thisBand] = [[] for j in range(len(regions[isFunded][thisParticle][thisBand].geoms))]

			# check for overlap with each polygon
			for j in range(0, len(regions[isFunded][thisParticle][thisBand].geoms)):
				if(regions[isFunded][thisParticle][thisBand].geoms[j].intersects(r)):
					contributingExps[isFunded][thisParticle][thisBand][j].append(experiment)
		
		# single connected region case
		else:
			if( not (thisParticle in contributingExps[isFunded])):
				contributingExps[isFunded][thisParticle] = {}
			if( not (thisBand in contributingExps[isFunded][thisParticle])):
				contributingExps[isFunded][thisParticle][thisBand] = []
			if(regions[isFunded][thisParticle][thisBand].intersects(r)):
				contributingExps[isFunded][thisParticle][thisBand].append(experiment)
		
		if(not isFunded): # if not funded also check for overlap with funded regions
		# multiple disconnected regions case

			if(regions[1][thisParticle][thisBand].geom_type == 'MultiPolygon'):
				if( not (thisParticle in contributingExps[1])):
					contributingExps[1][thisParticle] = {}
				if( not (thisBand in contributingExps[1][thisParticle])):
					contributingExps[1][thisParticle][thisBand] = [[] for j in range(len(regions[1][thisParticle][thisBand].geoms))]

				# check for overlap with each polygon
				for j in range(0, len(regions[1][thisParticle][thisBand].geoms)):
					if(regions[1][thisParticle][thisBand].geoms[j].intersects(r)):
						contributingExps[1][thisParticle][thisBand][j].append(experiment)
			
			# single connected region case
			else:
				if( not (thisParticle in contributingExps[1])):
					contributingExps[1][thisParticle] = {}
				if( not (thisBand in contributingExps[1][thisParticle])):
					contributingExps[1][thisParticle][thisBand] = []
				if(regions[1][thisParticle][thisBand].intersects(r)):
					contributingExps[1][thisParticle][thisBand].append(experiment)



print("plotting")
fig = plt.figure(figsize=[4.8*3, 4.8*2])
#gs = fig.add_gridspec(nrows=1, ncols=4, hspace=0, wspace=0, width_ratios=[9.02, 30.76715, 20, 12.55])
gs = gridspec.GridSpec(nrows=2, ncols=1, figure=fig, hspace=0.25)
gstop = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs[0], wspace=0)
gsbot = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs[1], wspace=0)
axs = [None]*2
axs[0] = gstop.subplots(sharex=False, sharey=False)
axs[1] = gsbot.subplots(sharex=False, sharey=True)
axs[0][0].axis('off')
axs[0][2].axis('off')
l = [None]*len(particles)
i = 0
for particle in particles:

	idIn = axIn[particle]
	idOut = axOut[particle]

	for band in regions[1][particle]:
		#patch = PolygonPatch(regions[0][particle], fc='white', ec=colors[particle], hatch=hatches[particle], alpha=0.25)
		#patch = PolygonPatch(regions[0][particle], fc=colors[particle], ec='dimgrey', hatch='/', alpha=0.25)
		#plt.gca().add_patch(patch)
		patch = PolygonPatch(regions[1][particle][band], fc=facecolors[particle], ec='dimgrey', alpha=0.25)
		#patch = PolygonPatch(regions[1][particle], fc='white', ec=colors[particle], hatch=hatches[particle], alpha=0.4)
		l[i] = axs[idOut][idIn].add_patch(patch)

	i += 1

for particle in particles:

	idIn = axIn[particle]
	idOut = axOut[particle]

	for band in regions[0][particle]:

		if(regions[0][particle][band].is_empty):
			continue

		patch = PolygonPatch(regions[0][particle][band], fc=(1,1,1,0.3), ec=edgecolors[particle], hatch=hatches[particle])#, alpha=0.8)
		#patch = PolygonPatch(regions[0][particle], fc=colors[particle], ec='dimgrey', hatch='/', alpha=0.25)
		axs[idOut][idIn].add_patch(patch)
		#patch = PolygonPatch(regions[1][particle], fc=colors[particle], ec='dimgrey', alpha=0.4)
		#plt.gca().add_patch(patch)


for i in range(0, 2):
	for j in range(0, 3):

		if(i == 0 and j != 1):
			continue

		axs[i][j].set_xscale('log')
		start, end = axs[i][j].get_xlim()
		lgstart = np.ceil(np.log10(start)/3)
		lgend = np.floor(np.log10(end)/3)
		lgticks = np.linspace(lgstart, lgend, 4)
		lgticks = 3*np.round(lgticks)
		axs[i][j].xaxis.set_ticks(10.**np.arange(3*lgstart, 3*lgend+1, 3))
		if(lgticks[0]-0.5 < np.log10(start)):
			#xmin = max(start/10, 10.**(lgticks[0]-1.0))#(lgticks[1]-lgticks[0])/3))
			xmin = 10.**(lgticks[0]-0.5)
		else:
			xmin = start*10.**0.3
		if(lgticks[-1]+0.5 >= np.log10(end)):
			xmax = 10.**(lgticks[-1]+0.5)
		else:
			xmax = end/10.**0.3
		axs[i][j].set_xlim(xmin, xmax)
		print(i, j, np.log10(xmin), np.log10(xmax), np.log10(xmax/xmin))
		widths[i,j] = np.log10(xmax/xmin)
		#axs[i].xaxis.set_ticks(10.**np.arange(lgstart, lgend, 3))
		#axs[i].xaxis.set_ticks(10.**lgticks)
		
		axs[i][j].set_ylim(ymin, ymax)

		#new_tick_locations = 10.**np.arange(-6,37,3)
		ax1 = axs[i][j]
		ax2 = ax1.twiny()
		ax2range = np.array(ax1.get_xlim())/hPlanck
		ax2.set_xlim(ax2range)
		ax2.set_xscale('log')
		start = ax2range[0]
		end = ax2range[1]
		lgstart = np.ceil(np.log10(start)/3)
		lgend = np.floor(np.log10(end)/3)
		#lgticks = np.linspace(lgstart, lgend, 4)
		lgticks = 3*np.floor(lgticks)
		#ax2.set_xticks(10.**lgticks[lgticks<30])
		ax2.set_xticks(10.**np.arange(3*lgstart, 3*lgend+1, 3))
		#if(j == 1):
			#ax1.set_xlabel('Energy (eV)')
			#ax2.set_xlabel('Frequency (Hz)')

# update width ratios so that x-axis aspect ratios are equal
topLR = (widths[1,:].sum()-widths[0,1])/2
widths[0,0] = topLR
widths[0,2] = topLR
gstop.set_width_ratios(widths[0,:])
gsbot.set_width_ratios(widths[1,:])
gs.update()


#ax2.set_xticks(new_tick_locations)
#ax2.set_xlabel(r"Frequency (Hz)")

#plt.xlim(1e-21, 1e22)
#plt.ylim(ymin, ymax)

#plt.xscale('log')
#plt.gca().xaxis.set_ticks(10.**np.arange(-21,22,3))
#plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%g'))
axs[0][1].yaxis.set_ticks(range(ymin, ymax+1, 5))
axs[1][0].yaxis.set_ticks(range(ymin, ymax+1, 5))

#fig.supxlabel('Energy (eV)')
#fig.supylabel('Year')
axs[0][1].set_ylabel('Year')
axs[1][0].set_ylabel('Year')
#fig.suptitle('Frequency (Hz)')
fig.text(0.512305, 0.07, 'Energy (eV)', ha='center')
fig.text(0.512305, 0.4825, 'Frequency (Hz)', ha='center')
fig.text(0.512305, 0.4975, 'Energy (eV)', ha='center')
fig.text(0.512305, 0.91, 'Frequency (Hz)', ha='center')

#new_tick_locations = 10.**np.arange(-6,37,3)
#iax1 = plt.gca()
#ax2 = ax1.twiny()
#ax2.set_xlim(np.array(ax1.get_xlim())/hPlanck)
#ax2.set_xscale('log')
#ax2.set_xticks(new_tick_locations)
#ax2.set_xlabel(r"Frequency (Hz)")


# fake line for legend
#l = [None]*len(particles)
lbls = [None]*len(particles)
i = 0
for particle in particles:
	#l[i] = plt.fill_between([1e0,1e0], 0, 0, facecolor=facecolors[particle], edgecolor='dimgrey', alpha=0.5)
	lbls[i] = particleLbls[particle]
	i += 1

fig.legend(handles=l,labels=lbls, loc=8, bbox_to_anchor=(0.5, 0.92, 0., 0.), numpoints=1, ncol=4, frameon=False, framealpha=0.)

# make labels and plot them
### no guarantees label placement will always work well ###
for isFunded in range(0, len(contributingExps)):

	# omit funded experiment labels
	#if(isFunded != 0):
	#	continue

	for thisParticle in contributingExps[isFunded]:

		idIn = axIn[thisParticle]
		idOut = axOut[thisParticle]

		for band in contributingExps[isFunded][thisParticle]:

			if(regions[isFunded][thisParticle][band].is_empty):
				continue
			
			# get label coordinates
			if(regions[isFunded][thisParticle][band].geom_type == 'MultiPolygon'):
				expLbl = [None]*len(regions[isFunded][thisParticle][band].geoms)
				x = [None]*len(regions[isFunded][thisParticle][band].geoms)
				y = [None]*len(regions[isFunded][thisParticle][band].geoms)
				for j in range(0, len(regions[isFunded][thisParticle][band].geoms)):
					thisx,thisy = regions[isFunded][thisParticle][band].geoms[j].exterior.coords.xy
					if(thisParticle == 'nu'):
						x[j] = np.asarray(thisx).max()/10.**0.1
						#if(band == 'PeV+'):
						#	x[j] = np.asarray(thisx).max()/10.**0.1
						#else:
						#	x[j] = np.asarray(thisx).min()*10.**0.1
					elif(thisParticle == 'cr' or (thisParticle=='gamma' and band == 'subMeV' and isFunded)):
						x[j] = np.asarray(thisx).min()*10.**0.1
					else:
						x[j] = np.asarray(thisx).max()/10.**0.1
					if(thisParticle == 'cr' or (thisParticle=='gamma' and band == 'subMeV' and isFunded)):
						y[j] = max(np.asarray(thisy).min(), ymin)
					else:
						y[j] = min(np.asarray(thisy).max(),ymax)
			else:
				expLbl = [None]*1
				x = [None]*1
				y = [None]*1
				thisx,thisy = regions[isFunded][thisParticle][band].exterior.coords.xy
				if(thisParticle == 'nu'):
					x[0] = np.asarray(thisx).max()/10.**0.1
					#if(band == 'PeV+'):
					#	x[0] = np.asarray(thisx).max()/10.**0.1
					#else:
					#	x[0] = np.asarray(thisx).min()*10.**0.1
				elif(thisParticle == 'cr'):
					if(isFunded and band == 'PeV+'):
						x[0] = np.asarray(thisx).min()*10.**0.1
					else:
						x[0] = np.asarray(thisx).max()/10.**0.1
				else:
					x[0] = np.asarray(thisx).max()/10.**0.1
				if(thisParticle == 'cr'):
					if(isFunded):
						y[0] = max(np.asarray(thisy).min(), ymin)
					else:
						y[0] = min(np.asarray(thisy).max(),ymax)
				else:
					y[0] = min(np.asarray(thisy).max(),ymax)
				#x[0], y[0] = regions[isFunded][thisParticle][band].representative_point().coords.xy
				#x[0] = x[0][0]
				#y[0] = y[0][0]
				"""
				while(not regions[isFunded][thisParticle][band].contains(sg.Point(x[0],y[0]))):
					if(thisParticle == 'gamma' and band == 'subMeV'):
						x[0] /= 10.**0.1
						y[0] -= 0.4
					else:
						x[0] /= 10.**0.1
						y[0] -= 0.01
					if(x[0] < 10.**-30):
						print("out of range")
						sys.exit(1)
				"""

			# make label
			polyList = contributingExps[isFunded][thisParticle][band]
			if(isinstance(polyList[0],list)): # check for multiple polygons
				for j in range(0, len(polyList)): # loop over polygons if there are multiple
					expLbl[j] = ""
					for k in range(0, len(polyList[j])):
						if(r'\n' in polyList[j][k]):
							splitString = polyList[j][k].split(r'\n')
							polyList[j][k] = ""
							for isub, sub in enumerate(splitString):
								polyList[j][k] += sub
								if(isub != len(splitString)-1):
									polyList[j][k] += "\n"
						if(k == 0):
							expLbl[j] += polyList[j][k]
						else:
							expLbl[j] += str(', %s' % polyList[j][k])

			else: # single polygon case
				expLbl[0] = ""
				for j in range(0, len(polyList)):
					if(r'\n' in polyList[j]):
						splitString = polyList[j].split(r'\n')
						polyList[j] = ""
						for isub, sub in enumerate(splitString):
							polyList[j] += sub
							if(isub != len(splitString)-1):
								polyList[j] += "\n"
					if(j == 0):
						expLbl[0] += polyList[j]
					else:
						expLbl[0] += str(', %s' % polyList[j])

			for j in range(0, len(expLbl)):
				
				""" FOR SOME REASON COORDINATES FROM TRANSFORM FUNCTIONS ARENT CORRECT, MAYBE SOMETHING ABOUT GRIDSPEC?
				if(thisParticle == 'cr'):
					_, ypixMax = axs[idOut][idIn].transData.transform((0,ymax))
				else:
					_, ypixMax = axs[idOut][idIn].transData.transform((0,y[j]))
				if(thisParticle=='nu' and band == 'PeV+'):
					_, ypixMin = axs[idOut][idIn].transData.transform((0,2026))
				else:
					_, ypixMin = axs[idOut][idIn].transData.transform((0,ymin))
				"""

				if(thisParticle=='nu'):
					#if(band == 'PeV+'):
					#	txt = axs[idOut][idIn].text(x[j],y[j], expLbl[j], color=textcolors[thisParticle], rotation=270, ha='right', va='top', fontsize=10, fontweight='bold', wrap=True)
					txt = axs[idOut][idIn].text(x[j],y[j], expLbl[j], color=textcolors[thisParticle], rotation=270, ha='right', va='top', fontsize=10, fontweight='bold', wrap=True)
					#else:
					#	txt = axs[idOut][idIn].text(x[j],y[j], expLbl[j], color=textcolors[thisParticle], rotation=270, ha='left', va='top', fontsize=10, fontweight='bold', wrap=True)
				elif(thisParticle=='cr'):
					if(not isFunded):
						txt = axs[idOut][idIn].text(x[j],y[j], expLbl[j], color=textcolors[thisParticle], rotation=270, ha='right', va='top', fontsize=10, fontweight='bold', wrap=True)
					elif(band == 'PeV+'):
						txt = axs[idOut][idIn].text(x[j],y[j], expLbl[j], color=textcolors[thisParticle], rotation=270, ha='left', va='bottom', fontsize=10, fontweight='bold', wrap=True)
					else:
						txt = axs[idOut][idIn].text(x[j],y[j], expLbl[j], color=textcolors[thisParticle], rotation=270, ha='right', va='bottom', fontsize=10, fontweight='bold', wrap=True)
				elif(thisParticle =='gamma' and band =='subMeV' and isFunded):
						txt = axs[idOut][idIn].text(x[j],y[j], expLbl[j], color=textcolors[thisParticle], rotation=270, ha='left', va='bottom', fontsize=10, fontweight='bold', wrap=True)
				else:
					txt = axs[idOut][idIn].text(x[j],y[j], expLbl[j], color=textcolors[thisParticle], rotation=270, ha='right', va='top', fontsize=10, fontweight='bold', wrap=True)

				if(thisParticle == 'nu'):
					if(band=='PeV+'):
						txtwidth=150
					else:
						txtwidth=250
				elif(thisParticle == 'cr'):
					if(band == 'PeV+'):
						txtwidth=120
					else:
						txtwidth=230
				elif(thisParticle == 'gamma'):
					if(band == 'MeV'):
						txtwidth = 75
					elif(band == 'subMeV'):
						if(isFunded):
							txtwidth = 150
						else:
							txtwidth = 75
					else:
						txtwidth = 230
				elif(thisParticle == 'gw' and isFunded):
					txtwidth = 180
				else:
					txtwidth = 230 
				txt._get_wrap_line_width = wrapwidth(txtwidth)

plt.savefig("gantt_chart_panelled2D.pdf", bbox_inches='tight')
