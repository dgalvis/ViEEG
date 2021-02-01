# It is generally advisable to run this

# Using anaconda, make the environment for python
#conda env create -f environmentVG_<os_name>.yml
# MATLAB must be installed

# This is only needed for mutual information methods (compiles the milca package)
./MI_initiator.sh

# The config.ini file is a subdirectory of ./net_results (NetRun.py and ModelRun.m already know that they need to look in net_results and analysis_results, so only enter subdirectory)
# Create Networks
ipython NetRun.py ./CC
ipython NetRun.py ./MI

# Run modelling
matlab -r "ModelRun('./CC')"
matlab -r "ModelRun('./MI')"

#All code was created by Daniel Galvis except where otherwise noted in README

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

