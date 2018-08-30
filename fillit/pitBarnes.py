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
from . import _fillBarnes
import numpy as np
import voropy as vpy

try: range = xrange
except: pass

class depressionFillingBarnes(object):
    """
    Perform pit filling algorithm based on Barnes 2014.

    Algorithm
    ---------
    Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and
    Watershed-Labeling Algorithm for Digital Elevation Models". Computers & Geosciences.
    Vol 62, Jan 2014, pp 117-127.
    """
    def __init__(self, Z=None, coords=None, cells=None, ngbIDs = None, ngbNb = None,
                        eps=0., sealevel=0., first=1, cartesian=True):
        """
        The class initialisation requires different variables depending of the type of grid used.

        Structured grid:
        -----------------
            + Z: 2D numpy array containing elevation values  (required)
            + eps: epsilon value in cases where a slope needs to be added to the flat depression  (optional)
            + sealevel: sea-level position, every nodes below sea-level will not be filled  (optional)
            + cartesian: needs to be set to .True. for structured grid  (optional)
            + first: if set to one will be used to compute the mesh parameters

        Untructured grid:
        -----------------
            + coords: 3D numpy array containing X,Y,Z coordinates  (required)
            + cells: 3D numpy array containing node IDs composing triangular cells (required)
            + ngbIDs: numpy array containing neighbors IDs for each node  (optional)
            + ngbNb: numpy array containing neighbors number for each node  (optional)
            + eps: epsilon value in cases where a slope needs to be added to the flat depression  (optional)
            + sealevel: sea-level position, every nodes below sea-level will not be filled  (optional)
            + cartesian: needs to be set to .False. for unstructured grid
            + first: if set to one will be used to compute the mesh parameters
        """

        self.cartesian = 0
        self.eps = eps

        if cartesian:
            self.cartesian = 1
            self.m = Z.shape[0]
            self.n = Z.shape[1]
        type = 0
        if self.eps > 0:
            type = 1

        # Regular grid initialisation
        if cartesian:
            _fillBarnes.fillinitialise(Z, type, first)
            self.Zin = Z
            self.seaIDs = np.where(Z<sealevel)
            return

        # Unstructured grid initialisation
        if first == 1:
            Tmesh = vpy.mesh_tri.MeshTri(coords, cells)
            boundary = Tmesh.get_boundary_vertices()
            edges_nodes = Tmesh.edges['nodes']
            coords = Tmesh.node_coords
            cells_nodes = Tmesh.cells['nodes']
            cells_edges = Tmesh.cells['edges']
            self.Zin = coords[:,2]
            self.m = len(self.Zin)
            self.seaIDs = np.where(self.Zin<sealevel)
            _fillBarnes.fillinitialise_unst_init(self.Zin, boundary, cells_nodes, cells_edges,
                                                        edges_nodes, type)
        elif first == 2:
            self.Zin = coords[:,2]
            self.m = len(self.Zin)
            self.seaIDs = np.where(self.Zin<sealevel)
            _fillBarnes.fillinitialise_unst_fast(self.Zin, boundary, ngbIDs, ngbNb, type)

        else:
            self.Zin = Z
            self.m = len(self.Zin)
            self.seaIDs = np.where(self.Zin<sealevel)
            _fillBarnes.fillinitialise_unst(self.Zin, type)

        return

    def getCellConnectivity(self):
        """
        Get cells connectivity for regular grid
        """

        nbcells = (self.m-1)*(self.n-1)

        return _fillBarnes.cellconnect(nbcells)

    def performPitFillingStruct(self):
        """
        Perform pit filling using Barnes' algorithm on structured grids!
        """

        if self.eps == 0.:
            fillZ = _fillBarnes.fillpit_struct(self.m,self.n)
        else:
            fillZ = _fillBarnes.fillpit_eps_struct(self.m,self.n,self.eps)

        fillZ[self.seaIDs] = self.Zin[self.seaIDs]

        return fillZ

    def performPitFillingUnstruct(self):
        """
        Perform pit filling using Barnes' algorithm on unstructured grids!
        """

        if self.eps == 0.:
            fillZ = _fillBarnes.fillpit_unstruct(self.m)
        else:
            fillZ = _fillBarnes.fillpit_eps_unstruct(self.m,self.eps)

        fillZ[self.seaIDs] = self.Zin[self.seaIDs]

        return fillZ
