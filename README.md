# octave-bg-sheath

This is the basic code, developed in GNUOctave, that I have used in my research in radial plasma sheath models. The model that I have developed is an extension of Allen-Boyd_Reynolds solution for the plasma sheath around a Langmuir probe, applied to a cylindrical Langmuir probe by Chen. The model was extended by Fernandez-Palop to include the ion temperature different from zero, but at that time the solution was only found for small enough ion temperature with respect to the electron temperature, which is the most frequent situation in cold plasmas. However, there was interest in obtaining a more general solution for several reasons:

*Electropositive plasmas have very few important applications, therefore negative ions should be included to model electronegative plasmas, which are much more important in industry.

*Under certain circumstances, the same electropositive Helium plasma behaves as predicted by ABR theory, but transitions to behave as predicted by the orbital complete theory of Berstein-Rabinowitz and Laframboise. In the search of a transitional theory between both ABR and BRL theories, an incomplete solution is of limited use. Therefore, a complete solution for any set of parameters, not only small ion temperature, is needed.

This code solves the system posed by Fernandez-Palop by imposing a very sensible condition, that is, that the physical magnitude cannot diverge, but this has some implications. In particular, when the ions have non zero temperature, the ion fluid behaves differently in sub-sonic and supersonic-regime, the system having a mathematical singularity when the ions reach the (local) speed of sound. This general statement was proposed by Valentini, but could not successfully be applied to this problem until now.

The code uses global variables. The initial reason, and the one that has precluded me from rewritting the code without globals, is that the Runge-Kutta method, or any other method for solving a differential equation system uses a function as an argument of the Runge-Kutta routine. The argument of that function has to be the independent variables of the system. In this case, the system has a set of constants, Ip, beta, alpha0, gamma (bg stands for beta-gamma, the most important constants) that could not be passed as parameters without defining a bigger set of independent variables some of which would be the constants and would not change. This global variable method makes the program simpler, so that it justifies its use.

This code is the basis of the following publications (I list the DOIs)
10.1063/1.4997844
10.1088/1361-6595/aaac58
10.1088/1361-6595/ab515e

And it is also important in these experimental publications
10.1088/1361-6587/ab3483
10.3390/app10175727

More routines have been used, most of them to sweep over values and get data for plotting. These have not been included.
