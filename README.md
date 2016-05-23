# PHYS227 HW 10

**Author:** Chinmai Raman

[![Build Status](https://travis-ci.org/chapman-phys227-2016s/hw-10-ChinmaiRaman.svg?branch=master)](https://travis-ci.org/chapman-phys227-2016s/hw-10-ChinmaiRaman)

**Due date:** 2016/04/28

## Specification

Complete the following exercises, placing your solutions into separate files. In each file, write the solution as a callable function, so that you can write suitable test functions that demonstrate correct output using the nose framework. GitHub will automatically run your tests on every commit, indicating any failures via the Travis framework with build status above.

1. Consider the coupled set of ODEs with initial conditions $u(0) = 1$, $v(0) = 0$:
   
   $$u'(t) = v(t)$$
   $$v'(t) = -u(t)$$
   
   Show that the solution $u(t) = \cos(t)$ and $v(t) = -\sin(t)$ satisfies these equations. Explain what $u(t)$ and $v(t)$ could correspond to in a physical system.
   
1. Using Euler's method, solve the above equations from $t=0$ until $t= 5(2\pi)$ using $N$ time steps per period (so $\Delta t = 2\pi / N$). What is the smallest integer value of $N$ that produces the qualitatively correct behavior of the solution over the solution interval? For 3 different choices of $N$, plot your simulated solution as a solid line on top of the exact solution as a dashed line for both $u(t)$ and $v(t)$. Be sure to annotate your plot with title, legend for the curves, and proper axes labels. Export the plots to png files in addition to loading them into your notebook.
   
1. Repeat using Heun's method.
   
1. Repeat using the 2nd-order Runge Kutta method.
   
1. Repeat using the 4th-order Runge Kutta method.

Finally, to cleanly present your work, create a Jupyter notebook ```hw10.ipynb``` that imports each of your python files as modules and demonstrates the functionality. The notebook should be formatted professionally using Markdown headings, including your name and a brief description of each homework section.

## Assessment

good assignment. Made me really learn the different methods of approximation.

## Honor Pledge

I pledge that all the work in this repository is my own with only the following exceptions:

* Content of starter files supplied by the instructor;
* Code borrowed from another source, documented with correct attribution in the code and summarized here.

Signed,

Chinmai Raman
