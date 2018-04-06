import sqlite3
import argparse
from cartopy import crs as ccrs
import matplotlib.pyplot as plt
import sys

def plotStatus(db):

    conn=sqlite3.connect(db)
    c = conn.cursor()
    c.execute("SELECT lon,lat,status FROM obs")
    data = c.fetchall()
    lons=list()
    lats=list()
    status=list()
    for row in data:
        lons+=[row[0]]
        lats+=[row[1]]
        status+=[row[2]]

    conn.close()

    # Plotting
    s=list()
    for o in range(0,len(status)):
        if status[o] <= 0:
            s += [40]
        else:
            s+=[60*int(status[o])]

    proj = ccrs.LambertConformal(central_longitude=15., central_latitude=63., standard_parallels=[63.])
    ax = plt.axes(projection=proj)
    ax.coastlines(resolution="10m")
    # Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r,
    # GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r,
    # Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r,
    # Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r,
    # Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r,
    # YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg,
    # brg_r, bwr, bwr_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r,
    # gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow,
    # gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray,
    # gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r,
    # ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spectral,
    # spectral_r, spring, spring_r, summer, summer_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
    plt.scatter(lons, lats, c=status, s=s, transform=ccrs.PlateCarree(),cmap="coolwarm")
    plt.colorbar()
    plt.show()

def plotProperty(db,prop):

    conn=sqlite3.connect(db)
    c = conn.cursor()
    c.execute("SELECT lon,lat,"+prop+" FROM obs where status <= 0")
    data = c.fetchall()
    lons=list()
    lats=list()
    values=list()
    for row in data:
        lons+=[row[0]]
        lats+=[row[1]]
        values+=[row[2]]

    conn.close()


    proj = ccrs.LambertConformal(central_longitude=15., central_latitude=63., standard_parallels=[63.])
    ax = plt.axes(projection=proj)
    ax.coastlines(resolution="10m")
    # Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r,
    # GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r,
    # Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r,
    # Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r,
    # Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r,
    # YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg,
    # brg_r, bwr, bwr_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r,
    # gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow,
    # gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray,
    # gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r,
    # ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spectral,
    # spectral_r, spring, spring_r, summer, summer_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
    plt.scatter(lons, lats, c=values, s=40, transform=ccrs.PlateCarree(),cmap="coolwarm",vmin=-5,vmax=5)
    plt.colorbar()
    plt.show()


def plotHistogram(db):

    conn=sqlite3.connect(db)
    c = conn.cursor()
    c.execute("SELECT lon,lat,fgdep,andep FROM obs where status <= 0 and fgdep > -30")
    data = c.fetchall()
    lons=list()
    lats=list()
    fgdeps=list()
    andeps=list()
    for row in data:
        lons+=[row[0]]
        lats+=[row[1]]
        fgdeps+=[row[2]]
        andeps+=[row[3]]

    conn.close()

    num_bins=200
    plt.subplot(2,1,1)
    plt.hist(fgdeps, num_bins, facecolor='blue', alpha=0.5)
    plt.subplot(2,1,2)
    plt.hist(andeps, num_bins, facecolor='red', alpha=0.5)
    #plt.colorbar()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Create parameter file for calibration methods (currently OI and override) in gridpp.i The program reads TITAN output and puts observed and CI values into the parameter file. If CI (rep column) is missing in the input, then a default value is used.')
    parser.add_argument('files', help='Netatmo input file', nargs="*")
    parser.add_argument('--db', help='Inputfile file', required=True)
    parser.add_argument('--type',help='Kind of plot',default="status")
    parser.add_argument('--prop', help='Poperty', default="")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    if args.type == "status":
        plotStatus(args.db)
    elif args.type == "prop":
        plotProperty(args.db,args.prop)
    elif args.type == "histogram":
        plotHistogram(args.db)
    else:
        raise NotImplementedError

if __name__ == "__main__":
    main()

