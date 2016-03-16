# A universal mathematical model of non-photochemical quenching

### Copyright (C) 2015-2016  Anna Matuszyńska, Oliver Ebenhöh

## Model
The mathematical model of the simplified photosynthetic electron transport chain to study non-photochemical quenching
and short-term light memory in higher plants. We provide here the code containing differential equations, reaction rates,
parameter space and initial conditions used for the simulations, together with the script to reproduce figures included
in the corresponding manuscript.

## Manuscript
Matuszynska et al.

## Experimental data
* DataBase/paperdata.db
* visualizeexperiment.py: visualize fluorescence and NPQ
* dataAnalysis.py: methods of retrieving information from teh database

## The core of the model
* cpfd.py: converts light intensity into light activation rate for specific specie
* lightProtocol.py: light protocols
* misc.py: pH, pHinv, colors
* npqModel.py: rates and equations
* npqResults.py: methods for analysis
* parameters.py: dictionary of the parameters
* simulate.py: integration (including steady state analysis)
* starter.py: contains all packages and instant. objects to run and play with the model from the console

##  Files reproducing manuscript figures
* reproducefigures.py: containd several functions named as fig1, fig4 or figs1, figs2
to reproduce figures from the manuscript with corresponding indexes

## Dependencies
All Python packages that are needed in order to run the model are listed in the file requirements.txt and can be installed by running

	pip install -r requirements.txt

It was noticed that under some systems LSODA integrator from the scipy.odeint package is not available.
This integrator is necessary for a correct computation. Please check if it is available on your system before running the model.

## How to start
To play with model from the console you can simply run the starter.py that imports all necessary modules
and creates all objects.

## Example
For the internal debbuging information we produce a number of output on the command line.
In order to avoid them, we suggest to suppress them by writing all log information to /dev/null.

To obtain similar graph as the Fig. 1 from the manuscript simply run:

    $ python
    >>> import reproducefigures
    >>> r = reproducefigures.ReproduceFigures()
    >>> r.fig1b()

To run create your own mutant, for instance without PTOX, you first need to import module with the parameters:

    >>> import parameters
    >>> p = parameters.NPQmodelParameterSet()
    >>> p.kPTOX = 0
    >>> m = npqModel.NPQModel(p, 2)

and than call the class than reproduce the figures with new, changed model as an argument

    >>> r = reproducefigures.ReproduceFigures(m)




