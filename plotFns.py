# Plot functions
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.pylab as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plotDef(rmax, p, xphm, yphm, phm):
    global colors
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rc("font", family="Helvetica Neue")
    plt.rc("xtick", labelsize="medium")
    linewidth = 1.2
    linestyle = "-"
    label_size = 15
    color = '#EE3984'
    lcolor = '#FFC708'

    fig, ax = plt.subplots(1,3,figsize=(12,4), gridspec_kw={'width_ratios': [1, 1, 1]}, dpi=150)
    ax[0].contourf(xphm, yphm, phm.c, cmap='Purples')
    thetPlt = np.linspace(0., np.pi, 200)
    ax[0].plot(p.a*np.cos(thetPlt), p.a*np.sin(thetPlt), c='white', zorder=1)

    ax[0].set_xlabel(r'$x$',usetex=True, fontsize=label_size)
    ax[0].set_ylabel(r'$y$',usetex=True,rotation=0, fontsize=label_size)
    ax[0].xaxis.set_label_coords(0.5,-0.1)
    ax[0].yaxis.set_label_coords(-0.15,0.4)
    ax[0].set_aspect(1)
    ax[0].set_ylim([-0.5,1.5])
    ax[0].set_xlim([-1.5,1.5])


    ax[1].set_xlabel(r'$t$',usetex=True, fontsize=label_size)
    ax[1].set_ylabel(r'$R(t)$',usetex=True,rotation=0, fontsize=label_size)
    ax[1].set_xlim([0., 10.])
    ax[1].set_ylim([0.,5.0])

    ax[1].xaxis.set_label_coords(0.5,-0.15)
    ax[1].yaxis.set_label_coords(-0.15,0.4)

    ax[1].spines["right"].set_visible(False)
    ax[1].spines["top"].set_visible(False)
    ax[1].spines["left"].set_position(("data", 0))
    ax[1].spines["bottom"].set_position(("data", 0))
    ax[1].set_aspect(1.3)

    ax[2] = plt.subplot(133, projection='polar')
    ax[2].set_thetalim(0, np.pi)
    ax[2].set_rlim(0, rmax)
    ax[2].set_xticks([])
    ax[2].set_xticklabels([])
    ax[2].set_yticks(np.linspace(0, rmax, 10))
    ax[2].set_yticklabels([])
    ax[2].tick_params("both", grid_alpha=0.50, grid_zorder=-10, grid_linewidth=0.5)
    

    # Theta ticks
    radius = ax[2].get_rmax()
    length = 0.025 * radius
    for i in range(180):
        angle = np.pi * i / 180
        ax[2].plot(
            [angle, angle],
            [radius, radius - length],
            linewidth=0.50,
            color="0.75",
            clip_on=False,
        )
    for i in range(0, 180, 5):
        angle = np.pi * i / 180
        ax[2].plot(
            [angle, angle],
            [radius, radius - 2 * length],
            linewidth=0.75,
            color="0.75",
            clip_on=False,
        )
    for i in range(0, 180, 15):
        angle = np.pi * i / 180
        ax[2].plot([angle, angle], [radius, 0], linewidth=0.5, color="0.75")
        ax[2].plot(
            [angle, angle],
            [radius + length, radius],
            zorder=500,
            linewidth=1.0,
            color="0.00",
            clip_on=False,
        )
        ax[2].text(
            angle,
            radius + 4 * length,
            "%d°" % i,
            zorder=500,
            rotation=i - 90,
            rotation_mode="anchor",
            va="top",
            ha="center",
            size="small",
            family="Helvetica Neue",
            color="black",
        )
    for i in range(0, 180, 90):
        angle = np.pi * i / 180
        ax[2].plot([angle, angle], [radius, 0], zorder=500, linewidth=1.00, color="0.0")


    # Radius ticks
    def polar_to_cartesian(theta, radius):
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        return np.array([x, y])


    def cartesian_to_polar(x, y):
        radius = np.sqrt(x ** 2 + y ** 2)
        theta = np.arctan2(y, x)
        return np.array([theta, radius])


    for i in np.linspace(0, rmax, rmax*2+1):
        P0 = 0, i
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, 0.75 * length]))
        ax[2].plot([P0[0], P1[0]], [P0[1], P1[1]], linewidth=0.50, color="0.75")

    for i in np.linspace(0, rmax, rmax+1):
        P0 = 0, i
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, +1.0 * length]))
        ax[2].plot([P0[0], P1[0]], [P0[1], P1[1]], zorder=500, linewidth=0.75, color="0.0")
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, -1.0 * length]))
        text = ax[2].text(
            P1[0],
            P1[1],
            "%d" % i,
            zorder=500,
            va="top",
            ha="center",
            size="x-small",
            family="Helvetica Neue",
            color="black",
        )
        text.set_path_effects(
            [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
        )
    return fig, ax


def plotDat(fig, ax, p, rxFull, ryFull, phxFull, phyFull, rwdFull, valFn, alph):
    # Plot data
    color = '#EE3984'
    lcolor = '#FFC708'
    colors = pl.cm.Spectral(np.linspace(0,1,p.nPtn))
    rcInd = int(p.nSc*p.l*p.nPtn/(np.pi*p.a))
    # global colors
    ax[0].scatter(rxFull, ryFull,
            s=15,
            color=color,
            # edgecolors='k',
            linewidth=0.5,
            alpha=alph/500,
            zorder=2)
    
    t = np.arange(np.size(rwdFull))*p.dt
    ax[1].plot(t, rwdFull,
            #    alpha=0.1,
               alpha=alph/10,
            #    color=colors[rcInd])
               color=lcolor)

    n = 2*p.nPtn + 1
    phiPlt = np.linspace(0, np.pi, n)
    # phiPlt[1::2] = phiPlt[0:-1:2]
    phiPlt[1:-2:2] = phiPlt[2:-1:2]
    phiPlt[-1] = phiPlt[0]
    R = np.zeros(n)

    R[1:-1:2] = valFn[:]
    R[:-1:2] = valFn[:]
    R[-1] = R[0]
    
    # R[:-2:2] = valFn[:]
    # R[1:-2:2] = valFn[1:]
    # R[-1] = R[0]
    # R[-2] = R[0]
    
    ax[2].fill(phiPlt, R, color="C1", zorder=150, alpha=alph/500)
    ax[2].plot(phiPlt, R, color="white", zorder=200, linewidth=1.5, alpha=alph/500)
    ax[2].plot(phiPlt, R, color="C1", zorder=250, linewidth=1, alpha=alph/500)
    
    plt.tight_layout()
    return fig, ax

def plotRadial(rmax, p, valFn):
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rc("font", family="Helvetica Neue")
    plt.rc("xtick", labelsize="medium")
    fig = plt.figure(figsize=(4, 4), dpi=150)
    ax = plt.subplot(1, 1, 1, projection="polar", frameon=True)
    ax.set_thetalim(0, np.pi)
    ax.set_rlim(0, rmax)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks(np.linspace(0, rmax, 10))
    ax.set_yticklabels([])
    ax.tick_params("both", grid_alpha=0.50, grid_zorder=-10, grid_linewidth=0.5)

    # Theta ticks
    radius = ax.get_rmax()
    length = 0.025 * radius
    for i in range(180):
        angle = np.pi * i / 180
        ax.plot(
            [angle, angle],
            [radius, radius - length],
            linewidth=0.50,
            color="0.75",
            clip_on=False,
        )
    for i in range(0, 180, 5):
        angle = np.pi * i / 180
        ax.plot(
            [angle, angle],
            [radius, radius - 2 * length],
            linewidth=0.75,
            color="0.75",
            clip_on=False,
        )
    for i in range(0, 180, 15):
        angle = np.pi * i / 180
        ax.plot([angle, angle], [radius, 0], linewidth=0.5, color="0.75")
        ax.plot(
            [angle, angle],
            [radius + length, radius],
            zorder=500,
            linewidth=1.0,
            color="0.00",
            clip_on=False,
        )
        ax.text(
            angle,
            radius + 4 * length,
            "%d°" % i,
            zorder=500,
            rotation=i - 90,
            rotation_mode="anchor",
            va="top",
            ha="center",
            size="small",
            family="Helvetica Neue",
            color="black",
        )
    for i in range(0, 180, 90):
        angle = np.pi * i / 180
        ax.plot([angle, angle], [radius, 0], zorder=500, linewidth=1.00, color="0.0")

    # Radius ticks
    def polar_to_cartesian(theta, radius):
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        return np.array([x, y])

    def cartesian_to_polar(x, y):
        radius = np.sqrt(x ** 2 + y ** 2)
        theta = np.arctan2(y, x)
        return np.array([theta, radius])


    for i in np.linspace(0, rmax, rmax*2+1):
        P0 = 0, i
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, 0.75 * length]))
        ax.plot([P0[0], P1[0]], [P0[1], P1[1]], linewidth=0.50, color="0.75")

    for i in np.linspace(0, rmax, rmax+1):
        P0 = 0, i
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, +1.0 * length]))
        ax.plot([P0[0], P1[0]], [P0[1], P1[1]], zorder=500, linewidth=0.75, color="0.0")
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, -1.0 * length]))
        text = ax.text(
            P1[0],
            P1[1],
            "%d" % i,
            zorder=500,
            va="top",
            ha="center",
            size="x-small",
            family="Helvetica Neue",
            color="black",
        )
        text.set_path_effects(
            [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
        )
    
    n = 2*p.nPtn + 1
    phiPlt = np.linspace(0, np.pi, n)
    # phiPlt[1::2] = phiPlt[0:-1:2]
    phiPlt[1:-2:2] = phiPlt[2:-1:2]
    phiPlt[-1] = phiPlt[0]
    R = np.zeros(n)

    R[1:-1:2] = valFn[:]
    R[:-1:2] = valFn[:]
    R[-1] = R[0]
    
    # R[:-2:2] = valFn[:]
    # R[1:-2:2] = valFn[1:]
    # R[-1] = R[0]
    # R[-2] = R[0]
    
    ax.fill(phiPlt, R, color="C1", zorder=150, alpha=0.1)
    ax.plot(phiPlt, R, color="white", zorder=200, linewidth=1.5, alpha=1.)
    ax.plot(phiPlt, R, color="C1", zorder=250, linewidth=1, alpha=1.)
    
    return fig, ax

def plotDefMul(rmax, p, xphm, yphm, phm):
    global colors
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rc("font", family="Helvetica Neue")
    plt.rc("xtick", labelsize="medium")
    linewidth = 1.2
    linestyle = "-"
    label_size = 15
    color = '#EE3984'
    lcolor = '#FFC708'

    fig, ax = plt.subplots(1,3,figsize=(12,4), gridspec_kw={'width_ratios': [1, 1, 1]}, dpi=150)
    ax[0].contourf(xphm, yphm, phm.c, cmap='Purples')
    thetPlt = np.linspace(0., np.pi, 200)
    ax[0].plot(p.a*np.cos(thetPlt), p.a*np.sin(thetPlt), c='white', zorder=1)

    ax[0].set_xlabel(r'$x$',usetex=True, fontsize=label_size)
    ax[0].set_ylabel(r'$y$',usetex=True,rotation=0, fontsize=label_size)
    ax[0].xaxis.set_label_coords(0.5,-0.1)
    ax[0].yaxis.set_label_coords(-0.15,0.4)
    ax[0].set_aspect(1)
    ax[0].set_ylim([-0.5,1.5])
    ax[0].set_xlim([-1.5,1.5])


    ax[1].set_xlabel(r'$t$',usetex=True, fontsize=label_size)
    ax[1].set_ylabel(r'$R(t)$',usetex=True,rotation=0, fontsize=label_size)
    ax[1].set_xlim([0., 10.])
    ax[1].set_ylim([0.,5.0])

    ax[1].xaxis.set_label_coords(0.5,-0.15)
    ax[1].yaxis.set_label_coords(-0.1,0.4)

    ax[1].spines["right"].set_visible(False)
    ax[1].spines["top"].set_visible(False)
    ax[1].spines["left"].set_position(("data", 0))
    ax[1].spines["bottom"].set_position(("data", 0))
    ax[1].set_aspect(1)
    
    ax[2] = plt.subplot(133, projection='polar')
    ax[2].set_thetalim(0, np.pi)
    ax[2].set_rlim(0, rmax)
    ax[2].set_xticks([])
    ax[2].set_xticklabels([])
    ax[2].set_yticks(np.linspace(0, rmax, 10))
    ax[2].set_yticklabels([])
    ax[2].tick_params("both", grid_alpha=0.50, grid_zorder=-10, grid_linewidth=0.5)


    # Theta ticks
    radius = ax[2].get_rmax()
    length = 0.025 * radius
    for i in range(180):
        angle = np.pi * i / 180
        ax[2].plot(
            [angle, angle],
            [radius, radius - length],
            linewidth=0.50,
            color="0.75",
            clip_on=False,
        )
    for i in range(0, 180, 5):
        angle = np.pi * i / 180
        ax[2].plot(
            [angle, angle],
            [radius, radius - 2 * length],
            linewidth=0.75,
            color="0.75",
            clip_on=False,
        )
    for i in range(0, 180, 15):
        angle = np.pi * i / 180
        ax[2].plot([angle, angle], [radius, 0], linewidth=0.5, color="0.75")
        ax[2].plot(
            [angle, angle],
            [radius + length, radius],
            zorder=500,
            linewidth=1.0,
            color="0.00",
            clip_on=False,
        )
        ax[2].text(
            angle,
            radius + 4 * length,
            "%d°" % i,
            zorder=500,
            rotation=i - 90,
            rotation_mode="anchor",
            va="top",
            ha="center",
            size="small",
            family="Helvetica Neue",
            color="black",
        )
    for i in range(0, 180, 90):
        angle = np.pi * i / 180
        ax[2].plot([angle, angle], [radius, 0], zorder=500, linewidth=1.00, color="0.0")


    # Radius ticks
    def polar_to_cartesian(theta, radius):
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        return np.array([x, y])


    def cartesian_to_polar(x, y):
        radius = np.sqrt(x ** 2 + y ** 2)
        theta = np.arctan2(y, x)
        return np.array([theta, radius])


    for i in np.linspace(0, rmax, rmax*2+1):
        P0 = 0, i
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, 0.75 * length]))
        ax[2].plot([P0[0], P1[0]], [P0[1], P1[1]], linewidth=0.50, color="0.75")

    for i in np.linspace(0, rmax, rmax+1):
        P0 = 0, i
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, +1.0 * length]))
        ax[2].plot([P0[0], P1[0]], [P0[1], P1[1]], zorder=500, linewidth=0.75, color="0.0")
        P1 = cartesian_to_polar(*(polar_to_cartesian(*P0) + [0, -1.0 * length]))
        text = ax[2].text(
            P1[0],
            P1[1],
            "%d" % i,
            zorder=500,
            va="top",
            ha="center",
            size="x-small",
            family="Helvetica Neue",
            color="black",
        )
        text.set_path_effects(
            [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
        )
    return fig, ax

def plotMultAgt(fig, ax, p, rxFull, ryFull, phxFull, phyFull, alph):
    # Plot data
    color = '#EE3984'
    lcolor = '#FFC708'
    acolor = '#7A7A7A'
    ax[0].scatter(rxFull, ryFull,
            # s=10,
            s=4,
            color=color,
            edgecolors='k',
            linewidth=0.5,
            alpha=alph,
            zorder=2)
    ax[0].quiver(rxFull, ryFull,
           phxFull, phyFull,
           scale_units='xy',
           color=acolor,
        #    scale=10.2)
           scale=30.2)
    return fig, ax

def plotMultFg(fig, ax, p, rxFull, ryFull, phxFull, phyFull, rwdFull, valFn, alph):
    # Plot data
    color = '#38BDF2'
    lcolor = '#FFC708'
    acolor = '#7A7A7A'
    ax[0].scatter(rxFull, ryFull,
            # s=10,
            s=4,
            color=color,
            edgecolors='k',
            linewidth=0.5,
            alpha=alph,
            zorder=2)
    ax[0].quiver(rxFull, ryFull,
           phxFull, phyFull,
           scale_units='xy',
           color=acolor,
           scale=30.2)
    
    t = np.arange(np.size(rwdFull))*p.dt
    ax[1].plot(t, rwdFull,
               alpha=alph/2,
               color=lcolor)

    n = 2*p.nPtn + 1
    phiPlt = np.linspace(0, np.pi, n)
    # phiPlt[1::2] = phiPlt[0:-1:2]
    phiPlt[1:-2:2] = phiPlt[2:-1:2]
    phiPlt[-1] = phiPlt[0]
    R = np.zeros(n)

    R[1:-1:2] = valFn[:]
    R[:-1:2] = valFn[:]
    R[-1] = R[0]
    
    # R[:-2:2] = valFn[:]
    # R[1:-2:2] = valFn[1:]
    # R[-1] = R[0]
    # R[-2] = R[0]
    
    ax[2].fill(phiPlt, R, color="C1", zorder=150, alpha=alph/500)
    ax[2].plot(phiPlt, R, color="white", zorder=200, linewidth=1.5, alpha=alph/500)
    ax[2].plot(phiPlt, R, color="C1", zorder=250, linewidth=1, alpha=alph/500)
    
    plt.tight_layout()
    # return fig, ax
    return fig, ax

def plotVidFn(p, xphm, yphm, phm, rxFull, ryFull, phxFull, phyFull, rwdFull):
        mpl.rcParams.update(mpl.rcParamsDefault)
        plt.rc("font", family="Helvetica Neue")
        plt.rc("xtick", labelsize="medium")
        linewidth = 1.2
        linestyle = "-"
        label_size = 15
        color = '#EE3984'
        lcolor = '#FFC708'
        acolor = '#7A7A7A'

        t = np.arange(np.size(rxFull))*p.dt
        # rwdFull = np.zeros(np.size(rxFull))
        fInd = 0
        for ind in np.arange(0,rxFull.shape[0],2):
                fig, ax = plt.subplots(1,2,figsize=(9,9), gridspec_kw={'width_ratios': [1, 1]}, dpi=150)
                ax[0].contourf(xphm, yphm, phm.c, cmap='Purples')
                thetPlt = np.linspace(0., np.pi, 200)
                ax[0].plot(p.a*np.cos(thetPlt), p.a*np.sin(thetPlt), c='white', zorder=1)
                ax[0].scatter(rxFull[ind], ryFull[ind],
                        s=20,
                        # s=4,
                        color=color,
                        edgecolors='k',
                        linewidth=0.5,
                        alpha=1,
                        zorder=2)
                ax[0].quiver(rxFull[ind], ryFull[ind],
                        phxFull[ind], phyFull[ind],
                        scale_units='xy',
                        color=acolor,
                        scale=10.2)
                        #  scale=30.2)
                ax[0].set_xlabel(r'$x$',usetex=True, fontsize=label_size)
                ax[0].set_ylabel(r'$y$',usetex=True,rotation=0, fontsize=label_size)
                ax[0].xaxis.set_label_coords(0.5,-0.05)
                ax[0].yaxis.set_label_coords(-0.1,0.5)
                ax[0].set_aspect(1)
                ax[0].set_ylim([-0.5,1.5])
                ax[0].set_xticks([])
                ax[0].set_yticks([])

                
                ax[1].plot(t[0:ind], rwdFull[0:ind], color = lcolor)
                ax[1].plot(t[ind], rwdFull[ind], 'o', color = lcolor)
                
                ax[1].set_xlabel(r'$t$',usetex=True, fontsize=label_size)
                ax[1].set_ylabel(r'$R(t)$',usetex=True,rotation=0, fontsize=label_size)
                ax[1].set_ylim([0., 5.0])
                ax[1].set_xlim([0., np.max(t)+2])

                ax[1].xaxis.set_label_coords(0.5,-0.15)
                ax[1].yaxis.set_label_coords(-0.15,0.5)

                ax[1].spines["right"].set_visible(False)
                ax[1].spines["top"].set_visible(False)
                ax[1].spines["left"].set_position(("data", 0))
                ax[1].spines["bottom"].set_position(("data", 0))
                ax[1].set_aspect(np.max(t)*0.7/5)
                plt.tight_layout()
                plt.savefig("Single"+str(fInd)+".png")
                plt.close()
                fInd += 1
        return

def plotMulVidFn(p, xphm, yphm, phm, rxFull, ryFull, phxFull, phyFull, rfgxFull, rfgyFull, phfgxFull, phfgyFull):
        mpl.rcParams.update(mpl.rcParamsDefault)
        plt.rc("font", family="Helvetica Neue")
        plt.rc("xtick", labelsize="medium")
        linewidth = 1.2
        linestyle = "-"
        label_size = 15
        color = '#EE3984'
        fgColor = '#38BDF2'
        lcolor = '#FFC708'
        acolor = '#7A7A7A'

        t = np.arange(np.size(rxFull))*p.dt
        rwdFull = np.zeros(np.size(rxFull))
        fInd = 0
        for ind in np.arange(0,rxFull.shape[0],2):
        # for ind in range(10):
                fig, ax = plt.subplots(1,2,figsize=(9,9), gridspec_kw={'width_ratios': [1, 1]}, dpi=150)
                ax[0].contourf(xphm, yphm, phm.c, cmap='Purples')
                thetPlt = np.linspace(0., np.pi, 200)
                ax[0].plot(p.a*np.cos(thetPlt), p.a*np.sin(thetPlt), c='white', zorder=1)
                ax[0].scatter(rxFull[ind], ryFull[ind],
                        # s=10,
                        s=4,
                        color=color,
                        edgecolors='k',
                        linewidth=0.5,
                        alpha=1,
                        zorder=2)
                ax[0].quiver(rxFull[ind], ryFull[ind],
                        phxFull[ind], phyFull[ind],
                        scale_units='xy',
                        color=acolor,
                        # scale=20.2)
                         scale=30.2)
                
                ax[0].scatter(rfgxFull[ind], rfgyFull[ind],
                        s=4,
                        color=fgColor,
                        edgecolors='k',
                        linewidth=0.5,
                        alpha=1,
                        zorder=2)
                ax[0].quiver(rfgxFull[ind], rfgyFull[ind],
                        phfgxFull[ind], phfgyFull[ind],
                        scale_units='xy',
                        color=acolor,
                         scale=30.2)
                
                ax[0].set_xlabel(r'$x$',usetex=True, fontsize=label_size)
                ax[0].set_ylabel(r'$y$',usetex=True,rotation=0, fontsize=label_size)
                ax[0].xaxis.set_label_coords(0.5,-0.05)
                ax[0].yaxis.set_label_coords(-0.1,0.5)
                ax[0].set_aspect(1)
                ax[0].set_ylim([-0.5,1.5])
                ax[0].set_xticks([])
                ax[0].set_yticks([])

                
                ax[1].plot(t[0:ind], rwdFull[0:ind], color = lcolor)
                ax[1].plot(t[ind], rwdFull[ind], 'o', color = lcolor)
                
                ax[1].set_xlabel(r'$t$',usetex=True, fontsize=label_size)
                ax[1].set_ylabel(r'$R(t)$',usetex=True,rotation=0, fontsize=label_size)
                ax[1].set_ylim([0., 1.0])
                ax[1].set_xlim([0., np.max(t)])

                ax[1].xaxis.set_label_coords(0.5,-0.15)
                ax[1].yaxis.set_label_coords(-0.15,0.5)

                ax[1].spines["right"].set_visible(False)
                ax[1].spines["top"].set_visible(False)
                ax[1].spines["left"].set_position(("data", 0))
                ax[1].spines["bottom"].set_position(("data", 0))
                ax[1].set_aspect(np.max(t)*0.7)
                plt.tight_layout()
                plt.savefig("Multi"+str(fInd)+".png")
                plt.close()
                fInd += 1
        return