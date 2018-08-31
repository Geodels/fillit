"""
Copyright 2018 Tristan Salles
fillIT is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.
fillIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License
along with fillIT.  If not, see <http://www.gnu.org/licenses/>.
"""

#!/usr/bin/python
# -*- coding: utf-8 -*-
from . import _fillZhou
import numpy as np
import meshplex as vpy

try: range = xrange
except: pass

class depressionFillingZhou(object):
    """
    Perform pit filling algorithm based on Zhou 2016.

    Algorithm
    ---------
    Zhou, Sun, Fu. "An efficient variant of the Priority-Flood algorithm for filling
    depressions in raster digital elevation models". Computers & Geosciences.
    Vol 90, Feb 2016, pp 87-96.
    """
    def __init__(self, Z=None, coords=None, cells=None, ngbIDs=None, ngbNb=None,
                        meshIDs=None, boundary=None, sealevel=0., cartesian=True,
                        extent=None, first=1):
        """
        The class initialisation requires different variables depending of the type of grid used.

        Structured grid:
        -----------------
            + Z: 2D numpy array containing elevation values  (required)
            + sealevel: sea-level position, every nodes below sea-level will not be filled  (optional)
            + cartesian: needs to be set to .True. for structured grid  (optional)
            + extent: specify if a partition is on the boundary of the domain [S,N,W,E] (optional)
            + first: if set to one will be used to compute the mesh parameters

        Untructured grid:
        -----------------
            + coords: 3D numpy array containing X,Y,Z coordinates  (required)
            + cells: 3D numpy array containing node IDs composing triangular cells (required)
            + ngbIDs: numpy array containing neighbors IDs for each node  (optional)
            + meshIDs: numpy array set to 1 for each node part of the pit filling mesh (optional)
            + ngbNb: numpy array containing neighbors number for each node  (optional)
            + boundary: numpy array containing local mesh boundary nodes  (optional)
            + sealevel: sea-level position, every nodes below sea-level will not be filled  (optional)
            + cartesian: needs to be set to .False. for unstructured grid
            + extent: numpy array set to 1 for each node on the boundary of the global domain  (optional)
            + first: if set to one will be used to compute the mesh parameters
        """

        # Regular grid initialisation
        if cartesian:
            self.m = Z.shape[0]
            self.n = Z.shape[1]
            self.Zin = Z
            if extent is None:
                extent = np.zeros((4),dtype=int)
            self.seaIDs = np.where(Z<sealevel)
            _fillZhou.fillinitialise(Z, extent, first)
            return

        if first > 0:
            self.Zin = coords[:,2]
        else:
            self.Zin = Z
        self.m = len(self.Zin)
        self.seaIDs = np.where(self.Zin<sealevel)
        if extent is None:
            extent = np.zeros((self.m),dtype=int)

        # Unstructured grid initialisation
        if first == 1:
            Tmesh = vpy.mesh_tri.MeshTri(coords, cells)
            boundary = Tmesh.get_boundary_vertices()
            edges_nodes = Tmesh.edges['nodes']
            coords = Tmesh.node_coords
            cells_nodes = Tmesh.cells['nodes']
            cells_edges = Tmesh.cells['edges']
            _fillZhou.fillinitialise_unst_init(coords, boundary, cells_nodes, cells_edges,
                                                        edges_nodes, extent)
        elif first == 2:
            _fillZhou.fillinitialise_unst_fast(coords, boundary, ngbIDs, ngbNb, meshIDs, extent)

        else:
            _fillZhou.fillinitialise_unst(self.Zin)

        return

    def performPitFillingStruct(self, simple=True):
        """
        Perform pit filling using Zhou's algorithm!
        """

        fillZ, pitlabs, watershed, graphnb = _fillZhou.fillpit_struct(self.m, self.n)

        fillZ[self.seaIDs] = self.Zin[self.seaIDs]

        if simple:
            return fillZ

        graph = _fillZhou.spillpts(graphnb)

        return fillZ, pitlabs, watershed, graph

    def performPitFillingUnstruct(self, simple=True):
        """
        Perform pit filling using Zhou's algorithm!
        """

        fillZ, pitlabs, watershed, graphnb = _fillZhou.fillpit_unstruct(self.m)

        fillZ[self.seaIDs] = self.Zin[self.seaIDs]

        if simple:
            return fillZ

        graph = _fillZhou.spillpts(graphnb)

        return fillZ, pitlabs, watershed, graph

    def getCellConnectivity(self):
        """
        Get cells connectivity for regular grid
        """

        nbcells = (self.m-1)*(self.n-1)

        return _fillZhou.cellconnect(nbcells)

    def combineRegTiles(self, updateFill, updateWatershed, localIDs):
        """
        Combines regular tiles when using parallel priority-flooding.
        """

        graph, nb = _fillZhou.combineregtiles(updateFill, updateWatershed, localIDs)

        return graph[:nb,:]

    def combineUnstructGrids(self, updateFill, updateWatershed, inids=None, outids=None):
        """
        Combines unstructured grids when using parallel priority-flooding.
        """

        if inids is not None:
            graph, nb = _fillZhou.combine_unstgrids(updateFill, updateWatershed, inids, outids)
        else:
            graph, nb = _fillZhou.combine_unstgrids_fast(len(outids), updateFill, updateWatershed)

        return graph[:nb,:]

    def fillGraph(self, cgraph, maxnghbs):
        """
        Perform pit filling over the spill-over graph.
        """

        return _fillZhou.graphfill(int(max(cgraph[:,1])+1), cgraph, maxnghbs)

    def getPitData(self, zi, zf, area, depID, totpit):
        """
        Extract pit information.
        """

        return _fillZhou.params_regpit(zi, zf, area, depID, totpit)

    def getPitData_unst(self, zi, zf, area, depID, totpit):
        """
        Extract pit information.
        """

        return _fillZhou.params_unstpit(zi, zf, area, depID, totpit)
