
import numpy as np
import pylab as pl
from matplotlib.ticker import FixedLocator, FixedFormatter
from scipy.spatial import cKDTree
from scipy.stats import nanmedian

import projmap
import globconnect
from hitta import WRY
import mycolor
import .figpref

import abstemp


#from njord import ghrsst
miv = np.ma.masked_invalid
an = lambda arr: arr[~np.isnan(arr) & ~np.isinf(arr)]
#oceanlabels = FixedFormatter(['Arctic', 'Atlantic', 'Indian',
#                              'Pacific', 'South'])
oceanlabels = FixedFormatter(['Arctic', 'Atl', 'Ind',
                              'Pacific', 'South'])


polyregs = np.array([
       5289, 5367, 5445, 5446, 5524, 5525, 5526, 5603, 5604, 5605, 5678,
       5679, 5680, 5749, 5750, 5818, 5819, 5820, 5885, 5886, 5887, 5888,
       5949, 5950, 5951, 5952, 6012, 6013, 6014, 6015, 6075, 6076, 6077,
       6078, 6136, 6137, 6138, 6139, 6196, 6197, 6198, 6199, 6200, 6257,
       6258, 6259, 6260, 6261, 6262, 6318, 6319, 6320, 6321, 6322, 6323,
       6376, 6377, 6378, 6379, 6380, 6381, 6436, 6437, 6438, 6439, 6440,
       6492, 6493, 6494, 6495, 6496, 6542, 6543, 6544, 6545, 6546, 6586,
       6587, 6588, 6589, 6590, 6628, 6629, 6630, 6631, 6632, 6666, 6667,
       6668, 6669, 6670, 6703, 6704, 6705, 6706, 6707, 6738, 6739, 6740,
       6741, 6742, 6743, 6773, 6774, 6775, 6776, 6777, 6778, 6809, 6810,
       6811, 6812, 6813, 6814, 6845, 6846, 6847, 6848, 6849, 6883, 6884,
       6885, 6886, 6887, 6923, 6924, 6925, 6926, 6927, 6972, 6973, 6974,
       6975, 6976, 7019, 7020, 7021, 7022, 7069, 7070, 7071, 7072, 7122,
       7123, 7124, 7171, 7172, 7173, 7216, 7218, 7262, 7264, 7265, 7316,
       7317, 7369])

def load():
    gl = globconnect.MinTimeMatrix('ecco025','full',regdi=8,regdj=8)
    gl.fix_mintmat()

    polyconn = an(gl.fixedmintmat[polyregs,:][:,polyregs])
    for r in polyregs:
        conn = lambda r: np.abs(an(gl.fixedmintmat[:,r]).shape)
        cur = conn(r)
        neg = conn(r-1) - cur
        pos = conn(r+1) - cur
        if (neg >= 50) | (pos >= 50):
            polyconn = gl.fixedmintmat[polyregs,r]
            if neg > pos:
                gl.fixedmintmat[:,r] = gl.fixedmintmat[:,r-1]
            else:
                gl.fixedmintmat[:,r] = gl.fixedmintmat[:,r+1]
            gl.fixedmintmat[polyregs,r] = polyconn


    gl.fixedmintmat[3537,:] = np.nan
    gl.fixedmintmat[:,3537] = np.nan
    gl.fixedmintmat[:,3604] = np.nan
    gl.fixedmintmat[3604,:] = np.nan
    gl.calc_dijkmat()

    gh = ghrsst.L4_K10()
    #gh.meansst = np.load('ghrsst_annual_sst.npz')['sst']
    gh.meansst = np.nanmax(np.load('ghrsst_month_sst.npz')['sst'], axis=0)

    gl.gcm.meansst = gl.gcm.reproject(gh, gh.meansst)
    gl.regsst = (np.bincount(gl.regmat.flat, gl.gcm.meansst.flat) /
                 np.bincount(gl.regmat.flat))
    #gl.sstmat = gl.regvec2grid(gl.regsst)
    return gl


def movedegree(degree=1):
    """Calculate the time it takes to move x degree"""
    gl = abstemp.open_mintmat_ds()
    mintmat = gl.dijkmintmat
    mintmat = mintmat.where(mintmat>=30, np.nan)
    gl.regdegvel = gl.regsst * np.nan
    for r in np.arange(gl.nreg):
        print r
        mask = gl.regsst >= (gl.regsst[r] + degree)
        try:
            gl.regdegvel[r] = an(mintmat[:,r][mask]).min()
        except ValueError:
            pass


def ssthist():

    jd1 = datestr2num('2010-01-01')
    jd2 = datestr2num('2013-12-31')
    jdvec = np.arange(jd1, jd2+1)
    hvec = np.linspace(5,40,100)
    hmat = np.zeros(jdvec.shape + hvec.shape)
    gh = ghrsst.L4_K10()
    for n,jd in enumerate(jdvec):
        try:
            gh.load(jd=jd)
            hmat[n,:-1],x = histogram(an(gh.sst), hvec)
            print "%i of %i" % (jd-jd1, jd2-jd1)
        except:
            pass
    return hmat




def tworegs_plot(gl):

    def reg(regid, dijk):
        tomat,frmat = gl.oneregmap(regid, dijk)
        idmat = gl.regmat * np.nan
        idmat[gl.regmat == regid] = 1
        xp,yp = gl.gcm.mp(gl.reglon[regid], gl.reglat[regid])
        return tomat,frmat,idmat,xp,yp

    def map(sp, mat, idmat, xp, yp, text):
        pl.subplot(4, 1, sp, axisbg="0.8")
        med = nanmedian(an(mat))
        gl.gcm.pcolor(mat, cmap=WRY(), rasterized=True, vmin=0, vmax=60)
        gl.gcm.mp.text(0, -75, "%s (%0.1f yrs)" % (text,med),
                       horizontalalignment="left")
        gl.gcm.pcolor(miv(idmat), rasterized=True)
        gl.gcm.mp.scatter(xp, yp, 80, facecolors='none', edgecolors='b', lw=1)

    def fig(F, dijk):
        pl.close(F)
        pl.figure(F, (4.98,8.82))
        pl.clf()
        pl.subplots_adjust(wspace=0,hspace=0)
        st = ['E', 'F', 'G', 'H'] if dijk else ['A', 'B', 'C', 'D']
        tomat,frmat,idmat,xp,yp = reg(800, dijk)
        map(1, tomat/365, idmat, xp, yp, '%s) Time to' % st[0])
        map(2, frmat/365, idmat, xp, yp, '%s) Time from' % st[1])
        tomat,frmat,idmat,xp,yp = reg(6900, dijk)
        map(3, tomat/365, idmat, xp, yp, '%s) Time to' % st[2])
        map(4, frmat/365, idmat, xp, yp, '%s) Time from' % st[3])
        mycolor.freecbar([0.12,0.07,0.78,0.010], [0,20,40,60], cmap=WRY(), size="large")
        pl.figtext(0.5, 0.04, 'Years', horizontalalignment="center", size="large")

    savesettings = {'dpi':300, 'bbox_inches':"tight"}
    fig(1, False)
    pl.savefig('figs/fig_3a_tofrom_examples_raw.png', **savesettings)
    pl.savefig('figs/fig_3a_tofrom_examples_raw.pdf', **savesettings)
    fig(2, True)
    pl.savefig('figs/fig_3b_tofrom_examples_dijk.png', **savesettings)
    pl.savefig('figs/fig_3b_tofrom_examples_dijk.pdf', **savesettings)



def mintmat(ax, gl, mat, tickregvec):
    ticks = [np.nonzero(tickregvec==p)[0].min() for p in [0,4,22,30,51]]
    labpos = (np.array(ticks + [len(tickregvec)])[1:] -
              np.array(ticks + [len(tickregvec)])[:-1]) / 2 + ticks
    pl.imshow(miv(mat ), cmap=WRY(), vmin=0, vmax=60)
    pl.xticks(ticks, [])
    pl.yticks(ticks, [])
    pl.grid(True, c="g", ls="-")
    pl.xlim(0, gl.nreg)
    pl.ylim(0, gl.nreg)
    ax.xaxis.set_minor_locator(FixedLocator(labpos))
    ax.xaxis.set_minor_formatter(oceanlabels)
    ax.yaxis.set_minor_locator(FixedLocator(labpos))
    ax.yaxis.set_minor_formatter(oceanlabels)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')

    pl.tick_params(which='minor', length=0,  bottom="off")
    pl.tick_params(which='major', direction="out", length=5, bottom="off",
                   colors="g",top="on",right="off")
    pl.setp(ax.yaxis.get_minorticklabels(), rotation=90, color="g")
    pl.setp(ax.xaxis.get_minorticklabels(), color="g")
    pl.ylabel("From region")
    pl.xlabel("To region")
    mycolor.freecbar([0.12,0.08,0.78,0.015], [0,20,40,60], cmap=WRY())
    pl.figtext(0.5, 0.04, 'Years', horizontalalignment="center", size="x-large")



def mintmat_plot_two_figs(gl):

    longvec   = longregvec(gl)
    longsort  = np.argsort(longvec)
    ticregvec = longvec[longsort]

    def fig(F):
        pl.close(F)
        pl.figure(F, (8,8))
        return pl.subplot(1,1,1,axisbg="0.5")

    longmintmat = gl.mintmat.astype(float)[:,longsort][longsort,:]
    longmintmat[longmintmat==0] = np.nan
    ax = fig(1)
    mintmat(ax, gl, longmintmat/365, ticregvec)
    pl.savefig('figs/fig_2a_mintime_raw.png', bbox_inches="tight", dpi=300)
    pl.savefig('figs/fig_2a_mintime_raw.pdf', bbox_inches="tight")

    longdijkmat = gl.dijkmintmat[:,longsort][longsort,:]
    longdijkmat[longdijkmat==0] = np.nan
    ax = fig(2)
    mintmat(ax, gl, longdijkmat/365, ticregvec)
    pl.savefig('figs/fig_2a_mintime_dijk.png', bbox_inches="tight", dpi=300)
    pl.savefig('figs/fig_2b_mintime_dijk.pdf', bbox_inches="tight")


def mintmat_plot_one_fig(gl):

    longvec   = longregvec(gl)
    longsort  = np.argsort(longvec)
    ticregvec = longvec[longsort]

    pl.close(1)
    pl.figure(1, (12, 6))
    pl.subplots_adjust(wspace=0.025)

    longmintmat = gl.mintmat.astype(float)[:,longsort][longsort,:]
    longmintmat[longmintmat==0] = np.nan
    ax = pl.subplot(1, 2, 1, axisbg="0.5")
    mintmat(ax, gl, longmintmat/365, ticregvec)
    pl.suptitle("Destination", size='x-large')
    ax.set_ylabel("Source", size='x-large')
    ax.set_xlabel("")
    ax.text(500,10116,"(A)", size="x-large")

    longdijkmat = gl.dijkmintmat[:,longsort][longsort,:]
    longdijkmat[longdijkmat==0] = np.nan
    ax = pl.subplot(1, 2, 2, axisbg="0.5")
    mintmat(ax, gl, longdijkmat/365, ticregvec)
    ax.set_xlabel("Dijkstra, to region")
    pl.setp( ax.yaxis.get_minorticklabels(), visible=False)
    pl.tick_params(which='major', length=0,  bottom="off")
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.text(500,10116,"(B)", size="x-large")

    pl.savefig('figs/fig_2_mintmat.png', bbox_inches="tight", dpi=300)
    pl.savefig('figs/fig_2_mintmat.pdf', bbox_inches="tight")






def longregvec(gl):

    fH = np.load('longhurst01.npz')
    llon = fH['llon']
    llat = fH['llat']
    larr = fH['array']

    kd = cKDTree(np.vstack((llon.flat,llat.flatten()+11.5)).T)
    dist,ij = kd.query(np.vstack((gl.reglon,gl.reglat)).T)
    longvec = np.zeros(len(ij))
    longvec[1:] = larr.flat[ij[1:]] + 1
    return longvec


def hists(gl):

    figpref.manuscript()
    pl.clf()

    yd,xd = np.histogram(an(gl.dijkmintmat/365), np.arange(0,100+1))
    xd = xd[:-1] + 0.5

    mintmat = gl.mintmat.copy()
    mintmat[mintmat==0] = np.nan
    yr,xr = np.histogram(an(mintmat/365), np.arange(0,100+1))
    xr = xr[:-1] + 0.5

    def plot():
        pl.step(xr, yr, where="mid", lw=2, c="b")
        pl.step(xd, yd, where="mid", lw=2, c="r")
        med = nanmedian(an(mintmat/365))
        pl.plot([med,med],pl.ylim(), ":", c="b")
        med = nanmedian(an(gl.dijkmintmat/365))
        pl.plot([med,med],pl.ylim(), ":", c="r")
        pl.legend(("Raw connectivty matrix", "Dijkstra connectivity matrix"))
        pl.ylabel("Number of connections")

    #pl.subplot(2,1,1)
    plot()
    pl.xlabel('Min Connectivity Time (Years)')
    pl.xlim(0,40)
    pl.savefig('figs/fig_4_mintime_hist.pdf', bbox_inches="tight")
    pl.savefig('figs/fig_4_mintime_hist.png', bbox_inches="tight")

    #pl.xlabel('')
    #pl.xlim(0,40)

    #pl.subplot(2,1,2)
    #plot()
    #pl.setp(pl.gca(), yscale="log")
    #pl.xlabel('Min Connectivity Time (Years)')
    #pl.savefig('figs/fig_4_mintime_hist_log.pdf', bbox_inches="tight")

def conn_region_map(gl):

    pl.close(1)
    pl.figure(1, (8.5,8))
    pl.subplots_adjust(wspace=0, hspace=0, top=0.925, bottom=0.125)

    pl.subplot(2,1,1)
    llon = np.load('data/Data_bror_community_X.py.npy')
    llat = np.load('data/Data_bror_community_Y.py.npy')
    clus = np.load('data/Data_bror_community_i.py.npy')
    mp = projmap.Projmap('glob')
    x,y = mp(llon,llat)
    mp.pcolormesh(x, y, clus,cmap=pl.cm.Paired, rasterized=True)
    mp.nice()
    gl.gcm.mp.text(0, -77, "A) Connectivity regions", horizontalalignment="left")


    pl.subplot(2,1,2)
    llon = np.load('data/Data_bror_centrality_X.py.npy')
    llat = np.load('data/Data_bror_centrality_Y.py.npy')
    clus = np.load('data/Data_bror_centrality_eig.py.npy')
    mp = projmap.Projmap('glob')
    x,y = mp(llon,llat)
    mp.pcolormesh(x, y, miv(clus), cmap=pl.cm.hot, rasterized=True)
    mp.nice()
    gl.gcm.mp.text(0, -77, "B) Eigenvector centrality", horizontalalignment="left")
    mycolor.freecbar([0.12,0.1,0.78,0.015], [0,0.25,0.5,0.75,1], cmap=pl.cm.hot)
    pl.figtext(0.5, 0.04, 'Centrality score',
               horizontalalignment="center", size="large")

    pl.savefig('figs/fig_5_connect_regions.pdf', bbox_inches="tight", dpi=300)
    pl.savefig('figs/fig_5_connect_regions.png', bbox_inches="tight", dpi=300)


def ecco_vel(gl):

    pl.close(1)
    pl.figure(1, (8,8))
    gl.gcm.load('uvel', jd=pl.datestr2num('2000-01-01'))
    gl.gcm.load('vvel', jd=pl.datestr2num('2000-01-01'))
    gl.gcm.pcolor(np.sqrt(gl.gcm.uvel**2 + gl.gcm.vvel**2)[0,:,:],
                  vmin=0, vmax=1, cmap=pl.cm.hot, rasterized=True)
    cb = pl.colorbar(orientation="horizontal", aspect=40, pad=0.020)
    cb.set_ticks([0,.25,0.5,0.75,1])
    cb.set_label("Water speed (m/s)")
    pl.savefig("figs/fig_1_ecco_waterspeed.pdf", bbox_inches="tight")
    pl.savefig("figs/fig_1_ecco_waterspeed.png", bbox_inches="tight", dpi=300)

def saturationplot(gl):

    mintmat = gl.mintmat.data
    tvec = np.unique(an(mintmat))
    satvec = [mintmat[mintmat==t].shape[0] for t in tvec]

    figpref.manuscript()
    pl.close('all')
    pl.plot(tvec/365,(np.cumsum(satvec)/123565456.)*100)
    pl.ylabel("Fraction of patches reached (%)")
    pl.xlabel("Time (Years)")
    pl.xlabel("Advection time (Years)")
    pl.savefig('figs/sup1_saturation.pdf',bbox_inches="tight")
    pl.savefig('figs/sup1_saturation.png',bbox_inches="tight")
