# Companion repository to *Computing submodules of points of general Drinfeld modules over finite fields*

This repository contains the SageMath implementation of the main algorithms of
the paper *Computing submodules of points of general Drinfeld modules over
finite fields*, by [Antoine
Leudière](https://cspages.ucalgary.ca/~antoine.leudiere1/) and [Renate
Scheidler](https://cspages.ucalgary.ca/~rscheidl/).

## Howto

All the functions are in the file [functions.sage](functions.sage). To use
them, the user can launch Sage and run
```
load('functions.sage')
```
provided that Sage is running in the same directory as the file
`functions.sage`.

However, we suggest to try our implementation by following the notebook
[notebook.ipynb](notebook.ipynb). There are two main ways to start the
notebook:
- Locally, the user may run the terminal command `sage -n jupyter` in the same
  directory as the `functions.sage` and `notebook.ipynb` files. The Jupyter
  interface will then launch in a web browser, and the user should click on the
  `notebook.ipynb` entry.
- Remotely, the user can try our notebook without any local installation of
  SageMath. To do so, simply click on
  [this link](https://mybinder.org/v2/gh/kryzar/research-drinfeld-submodules/HEAD).
  The user will be redirected to a *mybinder* webpage, which will deploy
  a running instance of SageMath in the brower. The initial loading will most
  likely take a few minutes.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kryzar/research-drinfeld-submodules/HEAD)

## Disclaimer

This is a proof-of-concept implementation, not ready for official SageMath
integration. If the algorithms prove valid, we plan to properly implement and
optimize them directly [in the main SageMath
codebase](https://github.com/sagemath/sage/pull/35026) (e.g. using the native
`DrinfeldModule_finite` class and the native implementation of Anderson
motives). See [this paper](https://dl.acm.org/doi/10.1145/3614408.3614417) and
[this
notebook](https://xavier.caruso.ovh/notebook/a_computational_approach_to_drinfeld_modules)
for the implementation of Drinfeld modules in SageMath.

The code is only tested on versions 10.7 and 10.8.
