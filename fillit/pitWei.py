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
from . import _fillWei
import numpy as np

try: range = xrange
except: pass

class depressionFillingWei(object):
    """
    Perform pit filling algorithm based on Wei 2018.

    Algorithm
    ---------

    Wei, Zhou, Fu. "Efficient Priority-Flood depression filling in raster digital
    elevation models", International Journal of Digital Earth. Jan 2018,
    DOI: 10.1080/17538947.2018.1429503
    """
    def __init__(self, Z, eps=0., cartesian=True):

        self.cartesian = 0
        self.eps = eps
        if cartesian:
            self.cartesian = 1
            self.m = Z.shape[0]
            self.n = Z.shape[1]
        else:
            print("Unable to  process Wei 2018 aklgorithm on irregular mesh")
            raise IOError('Irregular mesh unsupported for Wei 2018 algorithm...')

        type = 0
        _fillWei.fillinitialise(Z)

        return

    def performPitFilling(self):
        """
        Perform pit filling using Wei's algorithm!
        """

        return _fillWei.fillpit(self.m,self.n,self.eps)
