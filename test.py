import numpy as np

import getdist.plots as gplot

from getdist import plots
import pylab as plt
import numpy as np
from getdist import MCSamples
import os
import getdist.plots as gplot
import matplotlib.pyplot as plt
from scipy.misc import derivative, pade
from matplotlib import rc
plt.rcParams['text.usetex']=True
plt.rcParams.update({'font.size': 20})
dire=os.getcwd() # This will automatically set your current working directory
print(dire)



from getdist import MCSamples
from getdist import plots, MCSamples
import os
import getdist.plots as gplot

import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams['text.usetex']=True
plt.rcParams.update({'font.size': 20})
dire=os.getcwd() # This will automatically set your current working directory
print(dire)


import matplotlib as pl

pl.rcParams['figure.figsize'] = 14, 17
pl.rcParams['ytick.minor.visible'] =True
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['xtick.top'] = True
pl.rcParams['ytick.right'] = True
pl.rcParams['font.size'] = '20'
pl.rcParams['legend.fontsize'] = '10'
pl.rcParams['legend.borderaxespad'] = '1.9'
pl.rcParams['legend.numpoints'] = '1'

pl.rcParams['figure.titlesize'] = 'medium'
pl.rcParams['figure.titlesize'] = 'medium'
pl.rcParams['xtick.major.size'] = '10'
pl.rcParams['xtick.minor.size'] = '6'
pl.rcParams['xtick.major.width'] = '2'
pl.rcParams['xtick.minor.width'] = '1'
pl.rcParams['ytick.major.size'] = '10'
pl.rcParams['ytick.minor.size'] = '6'
pl.rcParams['ytick.major.width'] = '2'
pl.rcParams['ytick.minor.width'] = '1'
pl.rcParams['xtick.direction'] = 'in'
pl.rcParams['ytick.direction'] = 'in'
pl.rcParams['axes.labelpad'] = '10.0'
pl.rcParams['lines.dashed_pattern']=3.5, 1.0
#pl.rcParams['axes.formatter.limits']=-10,10
pl.rcParams['lines.dotted_pattern']= 1.0, 0.7

pl.rcParams['xtick.labelsize'] = '15'
pl.rcParams['ytick.labelsize'] = '15'
pl.rcParams['axes.labelsize'] = '30'
pl.rcParams['axes.labelsize'] = '30'

pl.rcParams['xtick.major.pad']='10'
pl.rcParams['xtick.minor.pad']='10'
#pl.rcParams['hatch.color'] = 'black'
pl.rc('axes', linewidth=2)

g1 = gplot.getSinglePlotter(chain_dir=r'./planck')
samples = g1.sampleAnalyser.samplesForRoot('base_plikHM_TTTEEE_lowl_lowE_lensing')
param = samples.getParams()
g1.settings.num_plot_contours = 2
g = plots.getSinglePlotter(width_inch=4, ratio=1)
P= param.sigma8*(param.omegam/0.3)**0.5
samples.addDerived(P, name=r'S80')#, label=r'c/h_rd')
k= param.H0
samples.addDerived(k, name=r'h')#, label=r'c/h_rd')
samples.updateBaseStatistics()
g.plot_2d([samples],[ r'h',r'S80'],filled=True,colors=['magenta'])


x=np.linspace(0,120,100)
plt.fill_between(x,0.70103627 ,0.77290325, color='darkgrey', alpha='0.5')
x=np.linspace(72.61,75.45,100)
plt.fill_between(x,0,5, color='turquoise', alpha='0.5')

names1 = ['Om0','rdrag','h','Or0','sigma8']
labels1 = [r'\Omega_{m0}',r'r_d',r'h',r'\Omega_{r0}',r'\sigma_8']
names2 = ['Om0','rdrag','w','h','Or0','sigma8']
labels2 = [r'\Omega_{m0}',r'r_d',r'w',r'h',r'\Omega_{r0}',r'\sigma_8']
names3 = ['Om0','rdrag','w0','wa','h','Or0','sigma8']
labels3 = [r'\Omega_{m0}',r'r_d',r'w_0',r'w_a',r'h',r'\Omega_{r0}',r'\sigma_8']
names4 = ['Om0','h','rdrag','P1','P2','Q1','Q2','sigma8']
labels4 =['Om0','h','rd','P1','P2','Q1','Q2','sigma8']
comb = [r'lcdm',r'wcdm',r'cpl',r'pade']
comb1=np.loadtxt('lcdm/lcdm_d+panth+Bao+fs8+Hz+5july.dat')
comb2=np.loadtxt('wcdm/wcdm_d+panth+Bao+Hz+fs8-5july.dat')
comb3=np.loadtxt('cpl/cpl_d+panth+Bao+Hz+fs8+5july.dat')
comb4=np.loadtxt('padesahintest/pade-d+Hz+panth+Bao+fs8+26june.dat')
print('Creating MCSamples for data.......')

c1         = MCSamples(samples=comb1,names = names1, labels = labels1)
c2         = MCSamples(samples=comb2,names = names2, labels = labels2)
c3         = MCSamples(samples=comb3,names = names3, labels = labels3)
c4         = MCSamples(samples=comb4,names = names4, labels = labels4)

param1 = c1.getParams()
P= param1.sigma8*(param1.Om0/0.3)**0.5
c1.addDerived(P, name=r'k')#, label=r'c/h_rd')
c1.updateBaseStatistics()
W=np.mean(P)
print(W)
o=c1.twoTailLimits(r'k',.68)
print(o)

print ([P.mean()])
print ([P.mean() -1.*P.std(),P.mean() +1.*P.std()])



"""
param2 = c2.getParams()
P= param2.sigma8*(param2.Om0/0.3)**0.5
c2.addDerived(P, name=r'k')#, label=r'c/h_rd')
c2.updateBaseStatistics()
W=np.mean(P)
print(W)
o=c2.twoTailLimits(r'k',.68)
print(o)

param3 = c3.getParams()
P= param3.sigma8*(param3.Om0/0.3)**0.5
c3.addDerived(P, name=r'k')#, label=r'c/h_rd')
c3.updateBaseStatistics()
W=np.mean(P)
print(W)
o=c3.twoTailLimits(r'k',.68)
print(o)

print('Pade')
param4 = c4.getParams()
P= param4.sigma8*(param4.Om0/0.3)**0.5
c4.addDerived(P, name=r'k')#, label=r'c/h_rd')
c4.updateBaseStatistics()


W=np.mean(P)
print(W)
o=c4.twoTailLimits(r'k',.68)
print(o)


param1 = c1.getParams()
k= param1.h*100.
c1.addDerived(k, name=r'H0')#, label=r'c/h_rd')
c1.updateBaseStatistics()

param2 = c2.getParams()
k= param2.h*100.
c2.addDerived(k, name=r'H0')#, label=r'c/h_rd')
c2.updateBaseStatistics()

param3 = c3.getParams()
k= param3.h*100.
c3.addDerived(k, name=r'H0')#, label=r'c/h_rd')
c3.updateBaseStatistics()

param4 = c4.getParams()
k= param4.h*100.
c4.addDerived(k, name=r'H0')#, label=r'c/h_rd')
c4.updateBaseStatistics()

c1.updateSettings({'contours': [0.683, 0.954,0.997]})
c2.updateSettings({'contours': [0.683, 0.954,0.997]})
c3.updateSettings({'contours': [0.683, 0.954,0.997]})
c4.updateSettings({'contours': [0.683, 0.954,0.997]})
g.settings.num_plot_contours = 2
g.plot_2d([c1,c2,c3,c4], r'H0', r'k',legend_labels=comb,filled=False,line_args=[{'ls':'-'},{'ls':'-.'},{'ls':':'},{'ls':'--'}])
g.add_legend([r'$Planck $',r'$ \Lambda CDM$',r'$wCDM$',r'$CPL$',r'$PADE$'], colored_text=True,align_right=True)


plt.ylabel(r'$ S_8$')
#plt.xlabel(r'$ \Omega_{m0} $')
plt.xlabel(r'$ H_{0}[kms^{-1} Mpc^{-1}] $')
g.export('H0s8kids.pdf') #LINE 20
"""
