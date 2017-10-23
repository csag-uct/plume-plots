#
# reads GCM data from historical and rcp files, merges them, calculates running median and CI
# plots time series of those with marking significant difference compared to a reference period 
# significance of differences calculated as pvalue of Mann-Whitney test
# CI calculated through bootstrapping differences in medians
# Piotr Wolski, May 2016
#
#

import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import pandas as pd
from matplotlib.colors import rgb2hex
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from datetime import timedelta, datetime
from netCDF4 import num2date, date2num

#from scipy.stats import mstats, mannwhitneyu, t, kendalltau
#import seaborn as sns
#from functools import partial
#from glob import glob
import sys
from statsmodels.distributions.empirical_distribution import ECDF

from arch.bootstrap import StationaryBootstrap, IIDBootstrap

#from scipy.stats import ttest_ind
#import statsmodels.stats.api as sm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.facecolor'] = 'eeeeee'
matplotlib.rcParams['axes.labelcolor'] = 'gray'
matplotlib.rcParams['axes.axisbelow'] = False
matplotlib.rcParams['axes.grid'] = False
matplotlib.rcParams['axes.linewidth'] = 0
matplotlib.rcParams['axes.titlesize'] = 10
matplotlib.rcParams['grid.color'] = 'white'
matplotlib.rcParams['grid.linestyle'] = '-'
matplotlib.rcParams['grid.linewidth'] = 2
matplotlib.rcParams['xtick.color'] = '666666'
matplotlib.rcParams['xtick.major.size'] = 0
matplotlib.rcParams['xtick.minor.size'] = 0
matplotlib.rcParams['ytick.color'] = '666666'
matplotlib.rcParams['ytick.major.size'] = 0
matplotlib.rcParams['ytick.minor.size'] = 0


#command line parameters
#datadir=sys.argv[1]
#varname=sys.argv[2]
#plotlabel=sys.argv[3]
#ylabel=sys.argv[4]
#regionname=sys.argv[5]
#p=int(sys.argv[6])
#statcode=sys.argv[7]

pdatadir="/home/lcoop/projects/WHO_hargrove/new_analysis/dscl_gridded/cmip5/year/timeseries/rcp45/pr/"

pdatadir="/terra/data/staging/local/users-a/lcoop/projects/Third_Nat_Com/data/cmip5/year/timeseries/historical/pr/"
fdatadir="/terra/data/staging/local/users-a/lcoop/projects/Third_Nat_Com/data/cmip5/year/timeseries/rcp85/pr/"
#pr_year_totals_MRI-CGCM3_historical_r1i1p1_1986-2005_SA.nc
varname="pr"
statcode="pr_year_totals"
plotlabel="Temporal evolution of total annual rainfall over southern Zambia \nin CMIP5 MME (rcp8.5)"
ylabel="pr [mm/year]"

#definitions of parameters
#gcms
gcms=['MIROC-ESM', 'CNRM-CM5', 'CanESM2', 'FGOALS-s2', 'BNU-ESM', 'MIROC5', 'MIROC-ESM-CHEM', 'MRI-CGCM3', 'bcc-csm1-1']
#gcms=['MIROC-ESM']
gcms=['inmcm4', 'MIROC-ESM', 'CNRM-CM5', 'CanESM2', 'FGOALS-s2', 'BNU-ESM', 'IPSL-CM5B-LR', 'MIROC5', 'HadGEM2-CC', 'MIROC-ESM-CHEM', 'MRI-CGCM3', 'MPI-ESM-LR', 'bcc-csm1-1', 'IPSL-CM5A-MR']

# reference period
fyear='1986'
lyear='2005'
nrefyears=20
#runing window - does not actually have to be the same as reference period, but haven't checked if that would work
nyears=20

#number of bootstraps
nboot=100
#if doing normal bootstrap (no autocorrelation present)
#bs=IIDBootstrap(np.arange(nyears))
#bsref=IIDBootstrap(np.arange(nrefyears))
#if doing block bootstrap (accounts for autocorrelation)
bs=StationaryBootstrap(3, np.arange(nyears))
bsref=StationaryBootstrap(3, np.arange(nrefyears))

#colors: [line_significant, shading_significant, line_nonsignificant, shading_nonsignificant]
cols=["red", "orange", "blue", "royalblue"]

#statcode="pr_year_totals"
#creating dataframe with input data
#p=5
#regionname="Winter rainfall region (Western Cape)"
#p=6
#regionname="Summer rainfall region (Gauteng)"

for gcm in gcms:
    pfile=pdatadir+statcode+"_"+gcm+"_historical_r1i1p1_1960-2005_SA.nc"
    ffile=fdatadir+statcode+"_"+gcm+"_rcp85_r1i1p1_2006-2099_SA.nc"
    print pfile
    ncdset=Dataset(pfile)
    pdat=ncdset.variables[varname][:]
    print  pdat.shape
    pdat=ncdset.variables[varname][:][:,(12,13,14),(6,7,8)]
    lat=ncdset.variables['lat'][:]
    lon=ncdset.variables['lon'][:]
    print  pdat.shape, lat, lon
    sys.exit()

    print pdat.shape
    pdat=np.mean(pdat,1)
    print pdat.shape
    pdat=np.mean(pdat,1)
    ptimes=ncdset.variables['time']
    pdates = num2date(ptimes[:],units=ptimes.units)
    pdat=np.ma.masked_invalid(pdat)

    ncdset=Dataset(ffile)
    fdat=ncdset.variables[varname][:]
    print fdat.shape
    fdat=np.mean(fdat,1)
    print fdat.shape
    fdat=np.mean(fdat,1)
    ftimes=ncdset.variables['time']
    fdates = num2date(ftimes[:],units=ftimes.units)
    fdat=np.ma.masked_invalid(fdat)

    alldat=np.append(pdat, fdat, 0)
    alldates=np.append(pdates,fdates,0)
    if gcm==gcms[0]:
       alldates0=np.copy(alldates)
#    gcms=[gcm+"_"+str(x) for x in pID]
#    print alldates
#    print alldates.shape, alldat.shape, len(gcmprov)
    datafr=pd.DataFrame(alldat, index=alldates0, columns=[gcm])
    if gcm==gcms[0]:
        alldata=datafr.copy()
    else:
        alldata=pd.concat([alldata,datafr], axis=1)

    ncdset.close()
#     print afile
#    ncdset=Dataset(afile)
    #this will pick up only one column from the nc file

#    lats=ncdset.variables['lat'][:]
#    lons=ncdset.variables['lon'][:]
#    latind=np.where(abs(lats-lat)==np.min(abs(lats-lat)))[0][0]
#    lonind=np.where(abs(lons-lon)==np.min(abs(lons-lon)))[0][0]
#    print latind, lonind
    #assumes data array is [lat,lon,time], adapt the line below if it isn't
#    print ncdset.variables[varname][:].shape
#    adata=ncdset.variables[varname][:][:,latind,lonind]

#    times=ncdset.variables['time']


#    dates = num2date(times[:],units=times.units)
#    adata=np.ma.masked_invalid(adata)

#    datafr=pd.DataFrame(adata, index=dates, columns=[gcm])
#    if gcm==gcms[0]:
#        alldata=datafr.copy()
#    else:
#        alldata=pd.concat([alldata,datafr], axis=1)
print alldata.shape


#processing and plotting

#crate figure object
fig=plt.figure(figsize=(10,5))
#and axis to plot to
pl = fig.add_subplot(1, 1, 1)

#prepare data
refdata=alldata[fyear:lyear]
fst=np.where(alldata.index==refdata.index[0])[0]
lst=np.where(alldata.index==refdata.index[-1])[0]
refindices=np.arange(fst,lst+1)

#reference period multimodel mean
allrefmean=refdata.mean().mean()
#these rescale all data so that each gcm has mean in reference period equal to multimodel mean
#allrefmean mean can be replaced by observed mean, if such is available
alldata=alldata-refdata.mean()+allrefmean
refdata=refdata-refdata.mean()+allrefmean

for gcm in gcms:
        gcmdata=alldata[gcm]
        refgcmdata=refdata[gcm]
        #yeah... magic happens here:
        #preparing arrays to store intermediate and final results
        N=len(gcmdata.values)
        ma=np.zeros(gcmdata.values.shape)
        ma[ma==0]=np.nan
        pval=np.copy(ma)
        ci_hi=np.copy(ma)
        ci_lo=np.copy(ma)
        subgcmdata=np.zeros([nyears,N])
        subgcmdata[subgcmdata==0]=np.nan
        indices=np.copy(subgcmdata)
        # running mean is calculated "by hand" and each window stored
        for i in range(nyears/2, N-nyears/2):
            subgcmdata[:,i]=gcmdata[i-nyears/2:i+nyears/2]
            indices[:,i]=np.arange(nyears)+i-nyears/2
        meandif=np.empty([nboot,N])
        b=0
        for bootind in bs.bootstrap(nboot):
            for refbootind in bsref.bootstrap(1):
                #this is to account for the overlap
                tempdata=subgcmdata[bootind[0][0],:]
                temprefdata=refgcmdata[refbootind[0][0]]
                #need to substitute if there is overlap
                for pos in np.arange(fst-nyears/2, lst+nyears/2)+1:
                    seldat=np.where(np.in1d(indices[:,pos],refindices))[0]
                    selref=np.where(np.in1d(refindices, indices[:,pos]))[0]
#                    temprefdata=np.copy(temprefdata0)
                    tempdata[seldat,pos]=temprefdata[selref]
                    #temprefdata[selref]=tempdata[seldat,pos]
                meandif[b,:]=np.nanmean(tempdata, 0)-np.nanmean(temprefdata)
                b=b+1
        #ci are calculated around the bootstrap mean. Otherwise, when bootstrap mean is biased, things might get ugly
        temp=np.percentile(meandif,50,0)
        ci_lo=np.percentile(meandif,5,0)-temp
        ci_hi=np.percentile(meandif,95,0)-temp
        for i in range(nyears/2, N-nyears/2):
            ma[i]=np.nanmean(subgcmdata[:,i])
            ecdf=ECDF(meandif[:,i])
            pv=ecdf(0)
            #to make sure it's one-tailed
            if pv>0.5:
               pv=1-pv
            #to make it two-tailed
            pv=pv*2
            pval[i]=pv
        if gcm==gcms[0]:
            MA=pd.DataFrame(ma, index=alldata.index, columns=[gcm])
            Pval=pd.DataFrame(pval, index=alldata.index, columns=[gcm])
            CI_lo=pd.DataFrame(ci_lo, index=alldata.index, columns=[gcm])
            CI_hi=pd.DataFrame(ci_hi, index=alldata.index, columns=[gcm])
        else:
            temp=pd.DataFrame(ma, index=alldata.index, columns=[gcm])            
            MA=pd.concat([MA, temp], axis=1)
            temp=pd.DataFrame(pval, index=alldata.index, columns=[gcm])            
            Pval=pd.concat([Pval, temp], axis=1)
            temp=pd.DataFrame(ci_lo, index=alldata.index, columns=[gcm])            
            CI_lo=pd.concat([CI_lo, temp], axis=1)
            temp=pd.DataFrame(ci_hi, index=alldata.index, columns=[gcm])            
            CI_hi=pd.concat([CI_hi, temp], axis=1)
#adjusting CI so they can be plotted
CI_hi=CI_hi+MA
CI_lo=CI_lo+MA

#plotting now
for gcm in gcms:
        emergeind=Pval[gcm]<0.05
        ind=np.invert(emergeind)
        indadj=np.copy(ind)
        for i in range(1,len(ind)-1):
            indadj[i]=np.sum(ind[i-1:i+2])
        indadj=np.invert(indadj)
        if np.sum(ind)>0:
            hi=CI_hi[gcm].copy()
            hi[ind]=np.nan
            lo=CI_lo[gcm].copy()
            lo[ind]=np.nan            
            pl.fill_between(alldata.index, lo, hi, color=cols[1], alpha=0.2)            
            hi=CI_hi[gcm].copy()
            hi[indadj]=np.nan
            lo=CI_lo[gcm].copy()
            lo[indadj]=np.nan            
            pl.fill_between(alldata.index, lo, hi, color=cols[3], alpha=0.2)            
#need to repeat to overplot ma on top of CI bands
for gcm in gcms:
        emergeind=Pval[gcm]<0.05
        ind=np.invert(emergeind)
        indadj=np.copy(ind)
        for i in range(1,len(ind)-1):
            indadj[i]=np.sum(ind[i-1:i+2])
        indadj=np.invert(indadj)
        if np.sum(ind)>0:
            mav=MA[gcm].copy()
            mav[ind]=np.nan
            pl.plot(mav, color=cols[0], linewidth=0.5, label='_nolegend_')
            mav=MA[gcm].copy()
            mav[indadj]=np.nan
            pl.plot(mav, color=cols[2], linewidth=0.5, label='_nolegend_')
pl.axhline(allrefmean, linewidth=1, color="black", label="reference period mean")
pl.set_title(plotlabel+"\n", fontsize=12, color="grey")
pl.set_xlabel("Year")
pl.set_ylabel(ylabel)
pl.spines['right'].set_visible(False)
pl.spines['top'].set_visible(False)
pl.xaxis.set_ticks_position('bottom')
pl.yaxis.set_ticks_position('left')        

line1 = mlines.Line2D([], [], color=cols[2], label='20-year moving average and its 95% confidence interval')
line2 = mlines.Line2D([], [], color=cols[0], label='as above, when difference wrt. reference period\'s mean significant at p=0.05')
line3 = mlines.Line2D([], [], color='black', linewidth=1, label='mean value in reference period (1986-2005)')

leg=plt.legend(handles=[line1, line2, line3], bbox_to_anchor=(0.655555, 0.3), ncol=1)
for text in leg.get_texts():
    text.set_color("grey")
    text.set_size(8)

fig.subplots_adjust(left=0.1, right=0.95, top=0.85, bottom=0.2, hspace=0.6, wspace=0.3)
plt.savefig("plume_"+statcode+".jpg")

