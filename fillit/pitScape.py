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
from . import _fillScape
import numpy as np
import meshplex as vpy

try: range = xrange
except: pass

class depressionFillingScape(object):
    """
    Perform pit filling algorithm for eSCAPE.
    """

    def __init__(self, Z=None, coords=None, cells=None, ngbIDs=None, ngbNb=None,
                 meshIDs=None, boundary=None, sealevel=0.,
                 extent=None, first=1):
        """
        The class initialisation requires different variables depending of the type of grid used.
        + coords: 3D numpy array containing X,Y,Z coordinates  (required)
        + cells: 3D numpy array containing node IDs composing triangular cells (required)
        + ngbIDs: numpy array containing neighbors IDs for each node  (optional)
        + meshIDs: numpy array set to 1 for each node part of the pit filling mesh (optional)
        + ngbNb: numpy array containing neighbors number for each node  (optional)
        + boundary: numpy array containing local mesh boundary nodes  (optional)
        + sealevel: sea-level position, every nodes below sea-level will not be filled  (optional)
        + extent: numpy array set to 1 for each node on the boundary of the global domain  (optional)
        + first: if set to one will be used to compute the mesh parameters
        """

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
            Tmesh.mark_boundary()
            ids = np.arange(0, len(Tmesh.node_coords), dtype=int)
            boundary = ids[Tmesh._is_boundary_node]
            coords = Tmesh.node_coords
            edges_nodes = Tmesh.edges['nodes']
            cells_nodes = Tmesh.cells['nodes']
            cells_edges = Tmesh.cells['edges']
            _fillScape.build_mesh(coords, boundary, cells_nodes, cells_edges,
                                           edges_nodes, extent)
        elif first == 2:
            _fillScape.escape_grid(coords, boundary, ngbIDs, ngbNb, meshIDs, extent)

        else:
            _fillScape.escape_grid_fast(self.Zin, 2)

        return

    def performPitFillingUnstruct(self, simple=True):
        """
        Perform pit filling!
        """

        fillZ, watershed, graphnb = _fillScape.fillpit(self.m)

        fillZ[self.seaIDs] = self.Zin[self.seaIDs]

        if simple:
            return fillZ

        graph = _fillScape.get_spillover_nodes(graphnb)

        return fillZ, watershed, graph

    def combineUnstructGrids(self, updateFill, updateWatershed, inids=None, outids=None):
        """
        Combines unstructured grids when using parallel priority-flooding.
        """

        if inids is not None:
            graph, nb = _fillScape.combine_edges(updateFill, updateWatershed, inids, outids)
        else:
            graph, nb = _fillScape.combine_edges_fast(len(outids), updateFill, updateWatershed)

        return graph[:nb,:]

    def fillGraph(self, cgraph, maxnghbs):
        """
        Perform pit filling over the spill-over graph.
        """

        elev,rank,nodes,spillID,order = _fillScape.global_graph_fill(int(max(cgraph[:,1])+2), cgraph, maxnghbs)
        order = order[::-1]
        keep = order > -1
        order = order[keep]
        
        graph = -np.ones((len(elev),4))
        graph[:,0] = elev
        graph[:,1] = nodes
        graph[:,2] = rank
        graph[:,3] = spillID

        return graph,order

    def getPitData_unst(self, zi, zf, area, depID, totpit):
        """
        Extract pit information.
        """

        return _fillScape.depression_info(zi, zf, area, depID, totpit)
