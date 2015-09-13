## Installation

To install the package simply type

	make install

The package requires *roxygen2*, *ggplot2*, *Matrix* and a couple of other libraries.

## Introduction

A first example:

	gp <- new.gp(0.5, kernel.exponential(1, 1))

This creates a Gaussian process with prior mean *0.5* and a squared exponential kernel. The likelihood model is a Gaussian with variance *1*. With

	gp <- new.gp(0.5, kernel.exponential(1, 1),
	      	     likelihood=new.likelihood("normal", 0.01))

the variance of the likelihood model is set to *0.01*. Assuming we have the following observations

	xp <- c(1, 2, 3)
	yp <- c(0.7, 0.7, 0.7)

we may compute the posterior Gaussian processs with

	gp <- posterior(gp, xp, yp)

The process can be summarized or plotted with

	summarize(gp, 1:100/20)

	plot(gp, 1:100/20)

### Link functions and other likelihood models

    TODO

### Heteroscedastic Gaussian process

    TODO
