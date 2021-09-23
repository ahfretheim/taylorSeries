# taylorSeries
A python library for computing and comparing taylor series expansions, approximations &amp; errors from taylor approximations. Can generate csv files with the comparative taylor constants, approximations or errors of an entire parametric family of functions specified by the user. More functionality likely to be added later. Note that this release is a Beta test.

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
