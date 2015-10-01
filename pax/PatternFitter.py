from __future__ import division
from collections import namedtuple
import json
import gzip

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin_powell
from scipy.ndimage.interpolation import zoom as image_zoom

from pax import utils
from pax.exceptions import CoordinateOutOfRangeException

# Named tuple for coordinate data storage
# Maybe works faster than dictionary... can always remove later
CoordinateData = namedtuple('CoordinateData', ('minimum', 'maximum', 'n_bin_edges', 'bin_spacing'))


class PatternFitter(object):

    def __init__(self, filename, zoom_factor=1, adjust_to_qe=None):
        """Initialize a pattern map file from filename.
        Format of the file is very similar to InterpolatingMap; a (gzip compressed) json containing:
            'coordinate_system' :   [['x', x_min, x_max, n_x], ['y',...
            'map' :                 [[valuex1y1, valuex1y2, ..], [valuex2y1, valuex2y2, ..], ...
            'name':                 'Nice file with maps',
            'description':          'Say what the maps are, who you are, your favorite food, etc',
            'timestamp':            unix epoch seconds timestamp
        Note x_min, x_max are bin *edges*, not bin centers.
        Likewise n_x is the number of x bin *edges*, so 1 larger than the number of bins in x.
        Sorry about that.

        zoom_factor is factor by which the spatial dimensions of the map will be upsampled

        adjust_to_qe is a list of same length as the number of pmts in the map;
            we'll adjust the patterns to account for these QEs, upweighing PMTs with higher QEs
            Obviously this should be None if map already includes QE effects (e.g. if it is data-derived)!
        """
        bla = gzip.open(utils.data_file_name(filename)).read()
        data = json.loads(bla.decode())

        self.data = np.array(data['map'])
        self.dimensions = len(data['coordinate_system'])    # Spatial dimensions (other one is sampling points)

        # Zoom the spatial map using linear interpolation, if desired
        if zoom_factor != 1:
            self.data = image_zoom(self.data, zoom=[zoom_factor] * self.dimensions + [1], order=1)

        # Adjust the expected patterns to the PMT's quantum efficiencies, if desired
        # No need to re-normalize: will be done in each gof computation anyway
        if adjust_to_qe is not None:
            self.data *= adjust_to_qe[[np.newaxis] * self.dimensions + [slice(None)]]

        # Store bin starts and distances for quick access, assuming uniform bin sizes
        self.coordinate_data = []
        for name, (start, stop, n_bin_edges) in data['coordinate_system']:
            n_bin_edges *= zoom_factor
            self.coordinate_data.append(CoordinateData(minimum=start,
                                                       maximum=stop,
                                                       n_bin_edges=n_bin_edges,
                                                       bin_spacing=(stop - start)/(n_bin_edges - 1)))

        # TODO: Technically we should zero the bins outside the tpc bounds again:
        # some LCE may have leaked into this region due to upsampling... but doesn't matter:
        # if it causes a bias, it will push some events who are already far outside the fiducial volume
        # even further out.

        self.n_points = self.data.shape[-1]
        self.default_point_selection = np.ones(self.n_points, dtype=np.bool)
        self.default_square_syst_errors = np.zeros(self.n_points)

    def expected_pattern(self, coordinates):
        """Returns expected pattern at coordinates -- NOT YET NORMALIZED!!!"""
        bes = self.get_bin_indices(coordinates)
        bes += [slice(None)]
        return self.data[bes].copy()

    def compute_gof(self, coordinates, areas_observed,
                    point_selection=None, square_syst_errors=None, statistic='chi2gamma'):
        """Compute goodness of fit at a single coordinate point
        :param areas_observed: arraylike of length n_points containing observed area at each point
        :param coordinates: arraylike of n_dimensions, coordinates to test
        :param point_selection: boolean array of length n_points, if False point will be excluded from statistic
        :param square_syst_errors: float array of length n_points, systematic error to use for each point
        :param statistic: 'chi2' or 'chi2gamma': goodness of fit statistic to use
        :return: value of goodness of fit statistic, or float('inf') if coordinates outside of range
        """
        return self._compute_gof_base(self.get_bin_indices(coordinates), areas_observed,
                                      point_selection, square_syst_errors, statistic)

    def compute_gof_grid(self, center_coordinates, grid_size, areas_observed,
                         point_selection=None, square_syst_errors=None, statistic='chi2gamma', plot=False):
        """Compute goodness of fit on a grid of points of length grid_size in each coordinate,
        centered at center_coordinates. All other parameters like compute_gof.
        Returns gof_grid, (bin number of lowest grid point in dimension 1, ...)
        :return:
        """
        bin_selection = []
        lowest_bins = []
        for dimension_i, x in enumerate(center_coordinates):
            cd = self.coordinate_data[dimension_i]
            if not cd.minimum <= x <= cd.maximum:
                raise CoordinateOutOfRangeException("%s is not in allowed range %s-%s" % (x, cd.minimum, cd.maximum))
            start = self._get_bin_index(max(x - grid_size / 2,
                                            self.coordinate_data[dimension_i].minimum),
                                        dimension_i)
            lowest_bins.append(start)
            stop = self._get_bin_index(min(x + grid_size / 2,
                                           self.coordinate_data[dimension_i].maximum),
                                       dimension_i)
            bin_selection.append(slice(start, stop + 1))        # Don't forget python's silly indexing here...

        gofs = self._compute_gof_base(bin_selection, areas_observed, point_selection, square_syst_errors, statistic)

        if plot:
            # Make the linspaces of coordinates along each dimension
            # Remember the grid indices are
            q = []
            for dimension_i, cd in enumerate(self.coordinate_data):
                dimstart = self._get_bin_center(bin_selection[dimension_i].start, dimension_i) - 0.5 * cd.bin_spacing
                # stop -1 for python silly indexing again...
                dimstop = self._get_bin_center(bin_selection[dimension_i].stop - 1, dimension_i) + 0.5 * cd.bin_spacing
                q.append(np.linspace(dimstart, dimstop, gofs.shape[dimension_i] + 1))
            q.append(gofs.T / np.nanmin(gofs))
            plt.pcolormesh(*q, vmin=1, vmax=4)
            plt.colorbar(label='Goodness of fit / minimum')

        return gofs, lowest_bins

    def get_bin_indices(self, coordinates):
        return [self._get_bin_index(x, dimension_i) for dimension_i, x in enumerate(coordinates)]

    def _get_bin_index(self, value, dimension_i):
        """Return bin index along dimension_i which contains value.
        Raises CoordinateOutOfRangeException if value out of range.
        TODO: check if this is faster than just using np.digitize on the bin list
        """
        cd = self.coordinate_data[dimension_i]
        if not cd.minimum <= value <= cd.maximum:
            raise CoordinateOutOfRangeException("%s is not in allowed range %s-%s" % (value, cd.minimum, cd.maximum))
        return int((value - cd.minimum) / cd.bin_spacing)

    def _get_bin_center(self, bin_i, dimension_i):
        cd = self.coordinate_data[dimension_i]
        return cd.minimum + cd.bin_spacing * (bin_i + 0.5)

    def _compute_gof_base(self, bin_selection, areas_observed, point_selection, square_syst_errors, statistic):
        """Compute goodness of fit statistic: see compute_gof
        bin_selection will be used to slice the spatial histogram.
        :return: gof with shape determined by bin_selection.
        """
        if point_selection is None:
            point_selection = self.default_point_selection
        if square_syst_errors is None:
            square_syst_errors = self.default_square_syst_errors
        square_syst_errors = square_syst_errors[point_selection]

        areas_observed = areas_observed.copy()[point_selection]
        total_observed = areas_observed.sum()
        fractions_expected = self.data[bin_selection + [point_selection]].copy()
        fractions_expected /= fractions_expected.sum(axis=-1)[..., np.newaxis]
        areas_expected = total_observed * fractions_expected

        # The actual goodness of fit computation is here...
        if statistic == 'chi2gamma':
            result = (areas_observed + np.clip(areas_observed, 0, 1) - areas_expected) ** 2
            result /= areas_expected + square_syst_errors + 1
        elif statistic == 'chi2':
            result = (areas_observed - areas_expected) ** 2
            result /= areas_expected + square_syst_errors
        else:
            raise ValueError('Pattern goodness of fit statistic %s not implemented!' % statistic)

        return np.sum(result, axis=-1)

    def minimize_gof_grid(self, center_coordinates, grid_size, areas_observed,
                          point_selection=None, square_syst_errors=None, statistic='chi2gamma', plot=False):
        """Return (spatial position which minimizes goodness of fit parameter, gof at that position)
        minimum is found by minimizing over a grid centered at center_coordinates
        and extending by grid_size in all dimensions
        All other parameters like compute_gof
        """
        gofs, lowest_bins = self.compute_gof_grid(center_coordinates, grid_size, areas_observed,
                                                  point_selection, square_syst_errors, statistic, plot)
        min_index = np.unravel_index(np.nanargmin(gofs), gofs.shape)
        # Convert bin index back to position
        result = []
        for dimension_i, i_of_minimum in enumerate(min_index):
            x = self._get_bin_center(lowest_bins[dimension_i] + i_of_minimum, dimension_i)
            result.append(x)

        if plot:
            plt.scatter(*[[r] for r in result], marker='*', s=20, color='orange', label='Grid minimum')

        return result, gofs[min_index]

    def minimize_gof_powell(self, start_coordinates, areas_observed,
                            point_selection=None, square_syst_errors=None, statistic='chi2gamma'):
        direc = None
        if self.dimensions == 2:
            # Hack to match old chi2gamma results
            s = lambda d: 1 if d < 0 else -1
            direc = np.array([[s(start_coordinates[0]), 0],
                              [0, s(start_coordinates[1])]])

        def safe_compute_gof(*args, **kwargs):
            try:
                return self.compute_gof(*args, **kwargs)
            except CoordinateOutOfRangeException:
                return float('inf')

        # Minimize chi_square_gamma function, fmin_powell is the call to the SciPy minimizer
        # It takes the function to minimize, starting position and several options
        # It returns the optimal values for the position (xopt) and function value (fopt)
        # A warnflag tells if the maximum number of iterations was exceeded
        #    warnflag 0, OK
        #    warnflag 1, maximum functions evaluations exceeded
        #    warnflag 2, maximum iterations exceeded
        rv = fmin_powell(safe_compute_gof,
                         start_coordinates, direc=direc,
                         args=(areas_observed, point_selection, square_syst_errors, statistic),
                         xtol=0.0001, ftol=0.0001,
                         maxiter=10, maxfun=None,
                         full_output=1, disp=0, retall=0)
        xopt, fopt, direc, iter, funcalls, warnflag = rv
        # On failure the minimizer seems to give np.array([float('inf')])
        if isinstance(fopt, np.ndarray):
            fopt = float('nan')
        return xopt, fopt
