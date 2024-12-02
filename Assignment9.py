import sys
sys.path.append('src')
import math

from fem import *
from collections import defaultdict
from elements.frame import FrameCreator
from elements.pointmass import PointmassCreator

# ----------------------------------------------------
#       mesh geometry: vertices, cells and nodesets
# ----------------------------------------------------

def createModel():
    vertices = [
        (0.0, 0.0),     # 0: Tree Base      # A
        (0.0, 3.0),     # 1                 # B
        (3.0, 4.0),     # 2: Light Bulb     # C
        (-3.0, 4.0),    # 3: Light Bulb     # D
        (-2.81, 3.41),  # 4: Light Bulb     # E
        (-2.27, 4.63),  # 5: Light Bulb     # F
        (-1.81, 3.24),  # 6: Light Bulb     # G
        (-1.0, 3.79),   # 7: Light Bulb     # H
        (1.77, 3.23),   # 8: Light Bulb     # I
        (1.0, 3.8),     # 9: Light Bulb     # J
        (2.78, 3.36),   # 10: Light Bulb    # K
        (2.15, 4.54),   # 11: Light Bulb    # L
        (-0.6, 6.59),   # 12: Light Bulb    # M
        (0.57, 6.59),   # 13: Light Bulb    # N
        (-1.8, 5.77),   # 14: Light Bulb    # O
        (1.78, 5.77),   # 15: Light Bulb    # P
        (-1.54, 5.08),  # 16: Light Bulb    # Q
        (1.36, 5.12),   # 17: Light Bulb    # R
        (-0.63, 5.8),   # 18: Light Bulb    # S
        (0.58, 5.84),   # 19: Light Bulb    # T
        (-0.81, 4.62),  # 20: Light Bulb    # U
        (0.77, 4.59),   # 21: Light Bulb    # V
        (0.0, 5.41),    # 22                # W
        (0.0, 4.0),     # 23                # Z
        (0.77, 3.25),   # 24                # A1
        (1.39, 3.46),   # 25                # B1
        (2.23, 3.74),   # 26                # C1
        (-1.01, 3.33),  # 27                # D1
        (-2.16, 3.72),  # 28                # E1
        (-0.45, 5.5),   # 29                # F1
        (0.25, 5.93)    # 30                # G1
    ]

    cells = [
        (0, 1),         # 0: Trunck
        (1, 23),        # 1: Trunck
        (23, 22),       # 2: Trunck
        (1, 24),        # 3: Limb
        (24, 25),       # 4: Limb
        (25, 26),       # 5: Limb
        (1, 27),        # 6: Limb
        (27, 28),       # 7: Limb
        (23, 20),       # 8: Limb
        (23, 21),       # 9: Limb
        (22, 29),       # 10: Limb
        (22, 30),       # 11: Limb
        (24, 9),        # 12: Branch
        (25, 8),        # 13: Branch
        (26, 10),       # 14: Branch
        (26, 2),        # 15: Branch
        (27, 7),        # 16: Branch
        (27, 6),        # 17: Branch
        (28, 4),        # 18: Branch
        (28, 3),        # 19: Branch
        (21, 11),       # 20: Branch
        (21, 17),       # 21: Branch
        (20, 16),       # 22: Branch
        (20, 5),        # 23: Branch
        (22, 15),       # 24: Branch
        (29, 18),       # 25: Branch
        (29, 14),       # 26: Branch
        (30, 19),       # 27: Branch
        (30, 13),       # 28: Branch
        (22, 12),       # 29: Branch
        (2,    ),       # 30: Light Bulb
        (3,    ),       # 31: Light Bulb
        (4,    ),       # 32: Light Bulb
        (5,    ),       # 33: Light Bulb
        (6,    ),       # 34: Light Bulb
        (7,    ),       # 35: Light Bulb
        (8,    ),       # 36: Light Bulb
        (9,    ),       # 37: Light Bulb
        (10,    ),      # 38: Light Bulb
        (11,    ),      # 39: Light Bulb
        (12,    ),      # 40: Light Bulb
        (13,    ),      # 41: Light Bulb
        (14,    ),      # 42: Light Bulb
        (15,    ),      # 43: Light Bulb
        (16,    ),      # 44: Light Bulb
        (17,    ),      # 45: Light Bulb
        (18,    ),      # 46: Light Bulb
        (19,    ),      # 47: Light Bulb
        (20,    ),      # 48: Light Bulb
        (21,    ),      # 49: Light Bulb
    ]

    vertexsets = {}
    vertexsets["Unions"] = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)

    cellsets = {}
    cellsets["Trunck"] = (0, 1, 2)
    cellsets["Limb"] = (3, 4, 5, 6, 7, 8, 9, 10, 11)
    cellsets["Branch"] = (12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29)
    cellsets["LightBulb"] = (30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49)

    # ----------------------------------------------------
    #       Create global data structures
    # ----------------------------------------------------
    theMesh = Mesh(vertices, cells, vertexsets, cellsets)
    theMesh.print()

    # ----------------------------------------------------
    # Model: elements types, constraints, and loading on the mesh
    # ----------------------------------------------------

    elmtTypes = {}
    g = 9.8
    E = 250e6 # Elastic Modulus Steel
    sigmae = 180e6 # Strength
    density = 7800 # Density

    R_Trunck = 30e-2
    A_Trunck = math.pi * (R_Trunck)**2
    I_Trunck = math.pi / 4 * R_Trunck**4
    W_Trunck = I_Trunck / R_Trunck

    R_Limb = 15e-2
    A_Limb = math.pi * (R_Limb)**2
    I_Limb = math.pi / 4 * R_Limb**4
    W_Limb = I_Limb / R_Limb

    R_Branch = 7.5e-2
    A_Branch = math.pi * (R_Branch)**2
    I_Branch = math.pi / 4 * R_Branch**4
    W_Branch = I_Branch / R_Branch

    elmtTypes["Trunck"] = FrameCreator({"E": E,
                                        "A": A_Trunck,
                                        "I": I_Trunck,
                                        "fx": 30,
                                        "fy": 0,
                                        "gravity": g,
                                        "density": density,
                                        "sigmae": sigmae,
                                        "w": W_Trunck})

    elmtTypes["Limb"] = FrameCreator({"E": E,
                                      "A": A_Limb,
                                      "I": I_Limb,
                                      "fx": 30,
                                      "fy": 0,
                                      "gravity": g,
                                      "density": density,
                                      "sigmae": sigmae,
                                      "w": W_Limb})

    elmtTypes["Branch"] = FrameCreator({"E": E,
                                        "A": A_Branch,
                                        "I": I_Branch,
                                        "fx": 30,
                                        "fy": 0,
                                        "gravity": g,
                                        "density": density,
                                        "sigmae": sigmae,
                                        "w": W_Branch})

    elmtTypes["LightBulb"] = PointmassCreator({"mass": 3000,
                                               "gravity": g})

    theModel = Model(theMesh, elmtTypes)
    theModel.addConstraint(vertexset = "Unions", dof = 0, value = 0.0)
    theModel.addConstraint(vertexset = "Unions", dof = 1, value = 0.0)
    theModel.addConstraint(vertexset = "Unions", dof = 2, value = 0.0)
    return theModel

def main():
    model = createModel()
    analysis = StaticLinearAnalysis(model)
    analysis.solve()
    model.printDetailed()


if __name__ == "__main__":
    main()