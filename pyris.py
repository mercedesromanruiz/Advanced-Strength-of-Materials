#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created on Thu Sep  1 16:49:17 2022
#
#@author: ignacio romero
#"""

import sys
sys.path.append('../../../src')


from examplebridge import createModel
from fem import StaticLinearAnalysis


# this line imports from examplebridge the structure
model = createModel()

# if we want to get information about the structure, uncomment the next line
# model.print()

# The next two lines take care of actually solving the structural problem
# The first line creates the analysis, the second one solves the equations
analysis = StaticLinearAnalysis(model)
analysis.solve()


# after the problem is solved, information about it can be printed on the screen
# This line prints on the screen information about the model and its solution
# uncomment it if you want detailed information
# model.printDetailed()
