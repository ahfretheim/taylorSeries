# -*- coding: utf-8 -*-
"""
taylorSeries library

A library of functions for analyzing infinitely differentiable functions 
according to their Taylor Series. Will use myself to analyze my sin^n - cos^n
function and then unleash upon the world. Also includes a recursive method for
calculating factorials.

REQUIRES: sympy & pandas (note: both are part of Anaconda)

IMPORTANT: DUE TO A QUIRK IN THE PYTHON PROGRAMMING LANGUAGE, YOU WILL NEED TO 
1) DEFINE YOUR SYMBOL IN YOUR LOCAL ENVIRONMENT/SCRIPT AND 2) LOCALLY IMPORT
BOTH SYMBOL AND SYMPY TO COMPLETE (1) AND WRITE YOUR EXPRESSIONS WHEN USING 
THIS LIBRARY.

Functions:
    
    factorial - Recursively computes the factorial of an integer. If a 
    non-integer number is given, truncates to integer. Has an optional input 
    of errorqty for negative numbers, as negative numbers do not have valid 
    factorials, which defaults to -1. If you instead wish to treat negative 
    numbers as positive numbers for your factorial, we reccomend that you use 
    Python's built-in abs function within the parentheses of the function to 
    convert.
    
    taylorConstants - Returns the first n taylor derivatives for the taylor
    series at a of a differentiable function as a list, with or without the
    factorial division applied. By default, this Taylor series is a Maclaurin
    series (a = 0), n defaults to 25 and factorialize is set to False.  Note 
    also that in this version of the functionality, only one variable can be
    defined, so the "symb" input should be a single variable (e.g. "n", "x"),
    and should also match the variable used in the expression.
    
    taylorApproximation - calculates the Nth order taylor approximation of your
    differentiable function.
    
    taylorApproximationList - recursively generates a list of ith-order Taylor
    approximations, in reverse order (e.g. highest order first). Primarily
    exists to support analyticSeries; I generally advise against direct calls
    to this function.
    
    errorfunction - Calculates the difference between the actual result at x,
    and a list of Taylor approximations. An optional boolean (defaulting to 
    False) enables all errors being given as absolute values, otherwise they
    are raw (default). Returns a list in the same order as the estimate list.
    
    analyticsToFile - A wrapper function for outputting analyticsSeries (below)
    to a csv file. Takes a filepath, and then all the inputs of analyticsSeries
    
    analyticSeries - Builds a pandas dataframe to compare functions varied by
    a parameter. Required inputs are the base expression, the x variable, and 
    the p variable to be varied for each column of the dataframe, in that
    order. Optional inputs include the order of the taylor polynomial (default
    is 25), the set of parameters to try (defaults to range 1 to 12), the A
    (defaults to Maclaurin series), a point to evaluate (defaults to 0.5), and 
    the mode (defaults to c, for constants). 
    
    A complete listing of modes is stated below:
        
        c - constants; the differentials of each order from 0 to 25
        a - approximation; the approximation of optional input X at each order
        f - factorialized constants; the constants divided by the proper
        factorial, as they would be in an actual Taylor expansion
        e - error; the error of the approximation at optional input X.
        ea - error absolute; error mode, but with absolute value applied
    
    Note that the optional X input is only relevant to the a, ea & e modes and 
    sees no use in c or f.
    
    TO ADD: function comparison methods & approximation comparison. Include
    analytic series function: takes TWO symbols, with one being varied for each
    column of analytic pandas dataframe.
"""

from sympy.core import symbols
from sympy import *
from pandas import DataFrame as df;

def analyticsToFile(filepath, expr, x, p, N = 25, pset = range(1, 12), A = 0, mode="c", X = 0.5):
    df = analyticSeries(expr, x, p, N, pset, A, mode, X)
    df.to_csv(filepath);
    return df; #currently returns dataframe for debug purposes, this will likely be removed after beta testing.

def analyticSeries(expr, x, p, N = 25, pset = range(1, 12), A = 0, mode="c", X = 0.5):
    L = {}
    for P in pset:
        E = expr.subs(p, P)
        if mode == "c":
            lst = taylorConstants(E, x, n = N, a = A)
            print(lst)
        elif mode == 'f':
            lst = taylorConstants(E, x, n = N, a = A, factorialize=True)
            print(lst)           
        elif mode == 'a':
            lst = taylorApproximationList(E, S = x, x = X, N = N, A = A)
            lst.reverse() #The recursive list generates the items in the order opposite of expected, so a mirroring method is used to correct this.
            print(lst)
        elif mode == 'e' or mode == 'ea':
            actual = E.subs(x, X)
            mst = taylorApproximationList(E, S = x, x = X, N = N, A = A)
            mst.reverse()
            if mode == 'e':
                abbie = False;
            elif mode == 'ea':
                abbie = True;
            else:
                abbie = False;
                print("WARNING: Logical leak detected in mode switch of analyticSeries method.")
            lst = errorfunction(mst, actual, ab=abbie)
            print(lst)
        else:
            print("Unrecognized mode. Returning trivial empty list.")
            lst = []
        L[str(E)] = lst
    return df.from_dict(L)

#recursive list building function for the a-mode functionality. Note that this list generates the items in reverse order:
def taylorApproximationList(E, S, x, N = 25, A = 0):
    n = int(N)
    ret = []
    print("Loop index " + str(n) + " commencing.")
    if n == 1:
        ret.append(taylorApproximation(E=E, S=S, x=x, N = n, A = A))
    elif n < 1:
        print("WARNING: Improper input.")
    else:
        ret.append(taylorApproximation(E=E, S=S, x=x, N = n, A = A))
        ret.extend(taylorApproximationList(E=E, S=S, x=x, N = N - 1, A = A))
    return ret

def taylorApproximation(E, S, x, N = 25, A = 0):
    ret = 0.0
    flist = taylorConstants(E, S, n = N, a = A, factorialize=True)
    if N > 1:
        r = range(0, N + 1)
    else:
        r = [0]
    for i in r:
        ret += flist[i]*(x - A)**i; #computes the ith polynomial term and adds to ret
    return ret

def taylorConstants(expr, symb, n = 25, a = 0, factorialize = False):
    ret = []
    ret.append(expr.subs(symb, 0))
    if n > 1:
        for i in range(1, n + 1):
            if factorialize:
                ret.append(float(diff(expr, symb, i).subs(symb, a))/factorial(i)); #determines ith derivative and then substitutes the value at a for the variable
            else:
                ret.append(diff(expr, symb, i).subs(symb, a)); #determines ith derivative and then substitutes the value at a for the variable;      
    elif n < 0:
        print("WARNING: Invalid quantity of taylor constants requested.")
    return ret
    

def factorial(n, errorqty = -1):
    N = int(n)
    if N == 1:
        return 1
    elif N == 0:
        return 1
    elif N < 0:
        return errorqty;
    else:
        return N*factorial(N-1)

def errorfunction(lst, actual, ab = False):
    ret = [];
    for L in lst:
        if ab:
            ret.append(abs(actual-L))
        else:
            ret.append(actual - L)
    return ret;