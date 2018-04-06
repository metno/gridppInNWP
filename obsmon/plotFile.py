
import netCDF4
import argparse
import numpy as np
from cartopy import crs as ccrs
import matplotlib.pyplot as plt
import sys

def readFile(fname,var):
    nc = netCDF4.Dataset(fname, 'r')
    print nc
    var_lons = nc["longitude"]
    var_lats = nc["latitude"]
    field = nc[var]
    print field.shape
    return nc,var_lons,var_lats,field

def closeFile(nc):
    nc.close()

def plotField(lons,lats,field,vmin,vmax,title):


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
    print field
    plt.title(title)
    plt.scatter(lons, lats, c=field, transform=ccrs.PlateCarree(),cmap="coolwarm",vmin=float(vmin),vmax=float(vmax),edgecolors='none')
    plt.colorbar()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('files', help='Netatmo input file', nargs="*")
    parser.add_argument('-f',dest="file",help='Inputfile file', required=True)
    parser.add_argument('-v',dest="var",help='variabel',required=True)
    parser.add_argument('--vmax', dest="vmax", default=5,help='vmax')
    parser.add_argument('--vmin', dest="vmin", default=-5, help='vmin')
    parser.add_argument('--title', dest="title", default="", help='Title')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    nc,lons,lats,field=readFile(args.file,args.var)
    print field.shape,lons.shape,lats.shape
    field=np.reshape(field,[lons.shape[0],lons.shape[1]])
    plotField(lons,lats,field,args.vmin,args.vmax,args.title)
    nc.close()

if __name__ == "__main__":
    main()

