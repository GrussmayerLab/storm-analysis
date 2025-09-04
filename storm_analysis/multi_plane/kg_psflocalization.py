#!/usr/bin/env python
"""
Giving a mapping file (from multi_plane.mapper), and a 
molecule list, generate molecule lists to use for the
PSF extraction step.

Hazen 05/17
"""

import numpy
import os
import pickle

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.sa_h5py as saH5Py


def psfCrossLocalizations(h5_filenames, frame = 0, aoi_size = 8, min_height = 0.0):

    # expand frames to length of filenames
    if isinstance(frame, int):
        t = frame
        frame = [t for _ in h5_filenames]

    # Load localizations & movie size.
    locs_mask = {}
    for i, filename in enumerate(h5_filenames):
        with saH5Py.SAH5Py(filename) as h5:
            locs = h5.getLocalizationsInFrame(frame[i])
            assert bool(locs), "No localizations found in frame " + str(frame[i])
            [movie_x, movie_y] = h5.getMovieInformation()[:2]

        # Remove localizations that are too dim.
        mask = (locs["height"] > min_height)

        locs_mask[i] = {}
        for elt in ["x", "y"]:
            locs_mask[i][elt] = locs[elt][mask]

        locs_mask[i]['fname'] = filename
    
    for k,v in locs_mask.items():
        # Remove localizations that are too close to each other.
        [xf, yf] = iaUtilsC.removeNeighbors(v["x"], v["y"], 2.0 * aoi_size)

        # Remove localizations that are too close to the edge or
        # outside of the image in any of the channels.
        #
        is_good = numpy.ones(xf.size, dtype = bool)
        for i in range(xf.size):
            
            # Check in Channel 0.
            if (xf[i] < aoi_size) or (xf[i] + aoi_size >= movie_x):
                is_good[i] = False
                continue
            
            if (yf[i] < aoi_size) or (yf[i] + aoi_size >= movie_y):
                is_good[i] = False
                continue

    #
    #
    # Save localizations for each channel.
    #
        gx = xf[is_good]
        gy = yf[is_good]

        basename = os.path.splitext(v['fname'])[0]
        saH5Py.saveLocalizations(basename + "_c1_psf.hdf5", {"x" : gx, "y" : gy})
    #
        # Print localizations that were kept.
        #
        print(gx.size, "localizations were kept out of", xf.size)
        for i in range(gx.size):
            print("ch{0}: {1:.2f} {2:.2f}".format(k, gx[i], gy[i]))

if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Determine localizations to use for PSF measurement.')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations file.")
    parser.add_argument('--map', dest='mapping', type=str, required=True,
                        help = "The name of the mapping file. This is the output of multi_plane.mapper.")
    parser.add_argument('--frame', dest='frame', type=int, required=False, default=0,
                        help = "The frame in .bin file to get the localizations from. The default is 0.")
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=8,
                        help = "The size of the area of interest around the bead in pixels. The default is 8.")
    parser.add_argument('--min_height', dest='min_height', type=float, required=False, default = 0.0,
                        help = "Minimum localization height.")

    args = parser.parse_args()
    
    psfLocalizations(args.mlist,
                     args.mapping,
                     frame = args.frame,
                     aoi_size = args.aoi_size,
                     min_height = args.min_height)
    
