import numpy as np
import math
import random


def equivalentStress(criterion, stress):
    """
    Compute the equivalent stress according to several yield
    criteria.

    Arguments:
    ----------
    criterion: a string with the name of the criterion. It can be "vonmises",
               "tresca"
    stress: a numpy array of dimensions 3x3

    Returns:
    ----------
    eq: the equivalent stress, a scalar
    """
    c = criterion.lower()

    if (c == "vonmises"):
        eq = (1/2 * ((stress[0, 0] - stress[1, 1])**2 + (stress[1, 1] - stress[2, 2])**2 + (stress[2, 2] - stress[0, 0])**2) + 3 * (stress[0, 1]**2 + stress[1, 2]**2 + stress[2, 0]**2))**(1/2)

    elif (c == "tresca"):
        eigenvalues = np.linalg.eigvals(stress)
        # print("Eigenvalues: ", eigenvalues)
        eq = np.max(eigenvalues) - np.min(eigenvalues)

    else:
        eq = 0.0

    return eq



def testStress():

    s = abs(random.random())
    stress = np.zeros((3,3))
    stress[0,0] = s
    print("1. Pure traction stress tensor:\n", stress)
    print("\n -- von Mises equivalent stress: ", equivalentStress("vonmises", stress))
    print(" -- von Mises error: ", s - equivalentStress("vonmises", stress))
    print(" -- Tresca equivalent stress", equivalentStress("tresca", stress))
    print(" -- Tresca error: ", s - equivalentStress("Tresca", stress))


    s = -abs(random.random())
    stress = np.zeros((3,3))
    stress[0,0] = s
    print("\n2. Pure compression stress tensor:\n", stress)
    print("\n -- von Mises equivalent stress: ", equivalentStress("vonmises", stress))
    print(" -- von Mises error: ", s + equivalentStress("vonmises", stress))
    print(" -- Tresca equivalent stress", equivalentStress("tresca", stress))
    print(" -- Tresca error: ", s + equivalentStress("tresca", stress))

    stress = np.zeros((3,3))
    s = abs(random.random())
    stress[0,1] = s
    stress[1,0] = s
    print("\n3. Pure shear stress tensor:\n", stress)
    print("\n -- von Mises equivalent stress", equivalentStress("vonmises", stress))
    print(" -- Tresca equivalent stress", equivalentStress("tresca", stress))
    print(" -- Tresca error: ", 2*s - equivalentStress("tresca", stress))

    m = np.random.rand(3, 3)
    stress = 0.5*(m+m.T)
    print("\n4. Random stress tensor:\n", stress)
    print("\n -- von Mises equivalent stress", equivalentStress("vonmises", stress))
    print(" -- Tresca equivalent stress", equivalentStress("tresca", stress))


if __name__ == "__main__":
    testStress()
