from sympy import *
from sympy import Rational as R
from IPython.display import display, Latex, HTML, Math
init_printing()
import numpy as np
import math
import pandas as pd


# this function will take the volume of tank a and b as well as the pipes
# and form the matrix
def tank(aTank, bTank, aPipe, bPipe, cPipe, dPipe, saltA, saltB):
    global A
    A = Matrix([[-R(aPipe + bPipe, aTank), R(bPipe, bTank)], [R(cPipe, aTank), -R(dPipe + bPipe, bTank)]])
    global y0
    y0 = Matrix([saltA, saltB])
    l = symbols('l')
    global l1, l2
    l1, l2 = solve(det(A - l * eye(np.shape(A)[0])))
    global v1
    v1 = (A - l1 * eye(np.shape(A)[0])).nullspace()[0]
    global v2
    v2 = (A - l2 * eye(np.shape(A)[0])).nullspace()[0]

    display("formed function")
    display(Math(r'A = ' + latex(A)))

    display("eigenvalues")
    display(A.eigenvals())
    display(Math(r'\lambda_1 =' + latex(l1) + r'\approx' + latex(round(l1, 2))))
    display(Math(r'\lambda_2 =' + latex(l2) + r'\approx' + latex(round(l2, 2))))

    display("eigenvectors")
    display(A.eigenvects())
    display(Math(r'v_1 = ' + latex(v1) + r'=' + latex(v1.evalf(4))))
    display(Math(r'v_2 = ' + latex(v2) + r'=' + latex(v2.evalf(4))))

    display("y0")
    display(y0)

    global c1
    c1 = Matrix.hstack(v1, v2, y0).rref()[0][0, -1]
    global c2
    c2 = Matrix.hstack(v1, v2, y0).rref()[0][1, -1]
    display("c1 and  c2")
    display(c1, c2)

    return

    # this fuction takes the time and retruns the vector with the salt
    # after the time specified in the params. First entry in the output
    #  is for tank A, second is for B. The output might seem weird,
    #  therefore is good to call evalf() on it.


def saltInTankAtTime(time):
    display(c1 * v1 * exp(time * l1) + c2 * v2 * exp(time * l2))
    return (c1 * v1 * exp(time * l1) + c2 * v2 * exp(time * l2))


# this function takes the max amout
# of time and the fraction of salt that will
# remain from the initial amount, It will output the
# specific time when the target is met, and the amount
# of salt reached
def findTimeWhenFractionOfInitial(tMax, fraction):
    C = v1.row_join(v2).row_join(y0).rref()[0][:, -1]

    for t in range(0, tMax):
        y1 = float(C[0] * v1[0] * math.e ** (l1 * t) + C[1] * v2[0] * math.e ** (l2 * t))
        if y1 <= fraction:
            display(t)
            display(y1)
            break
    return


# this fuctuon takes the min time, the max time
# and an amount that represents the salt, and will
# output the time when the two tanks have a difference
# less than the specified amount
def findTimeWhenDifferenceIsLessThan(tMin, tMax, difference):
    for x in range(tMin, tMax):
        f1 = c1 * v1[0] * math.e ** (l1 * x) + c2 * v2[0] * math.e ** (l2 * x)
        f2 = c1 * v1[1] * math.e ** (l1 * x) + c2 * v2[1] * math.e ** (l2 * x)
        if f1 - f2 <= difference:
            print(x)
            break


if __name__ == '__main__':
    tank(200,100,5,8,3,5,100,0)