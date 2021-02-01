#All code was created by Daniel Galvis except where otherwise noted (MEGmat/BNI_find.m, MEGmat/BNI_single.m MEGmat/NI_model.m, MEGmat/theta_model.m, MEGmat/bonf_holm.m, MEGpy/MutualInformationICA/*, and MEGpy/connectivity.py->out=MIhigherdim(d1,d2))
#
#Copyright (C) 2019 Daniel Galvis
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

# importing MEGpy imports MEGpy.preprocess and MEGpy.connectivity and MEGpy.connectivity_v2.py
# Be careful since connectivity and connectivity_v2 have similarly named files!
from MEGpy import preprocess
from MEGpy import connectivity
__all__ = ['preprocess', 'connectivity']

