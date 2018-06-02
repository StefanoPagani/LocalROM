# LocalROM

**MATLAB(R) scripts for the numerical approximation of the Fitzhugh-Nagumo membrane model using local reduced-order models.**

LocalROM has been developed at [MOX-Politecnico di Milano](https://mox.polimi.it) under the support of the ERC Advanced Grant [iHEART project](http://iheart.polimi.it). It includes straightforward implementations of the local reduced-order models presented in the submitted article:
>[**[PMQ18] S. Pagani, A. Manzoni, A. Quarteroni. Numerical approximation of parametrized problems in cardiac electrophysiology by a local reduced basis method**, submitted to CMAME, 2018.]


Download and Installation
-------

To install the library, extract the ZIP file or clone the git repository.

Run the script by running the setup file
```Matlab
setPath
```

Examples
-------

* **FOM**:
numerical approximation of the parametrized Fitzhugh-Nagumo membrane model for an instance of the parameter epsilon.

![Full-order model approximation](/Figures/FOMexample.png)

* **GlobalROMConvergence**:
mean relative error convergence analysis with respect to the number of POD basis functions for the global reduced-order model.

* **GlobalROMHyperred**:
mean relative error convergence analysis with respect to the number of POD basis functions for the global hyperreduced model.

* **LocalROM**:
mean relative error convergence analysis with respect to the number of POD basis functions for the local reduced-order models (time-, parameter- and state-based).

![Local ROM error convergence](/Figures/LocalROMConvergence.png)

* **LocalROMHyperred**:
mean relative error convergence analysis with respect to the number of POD basis functions for the local hyperreduced models (time-, parameter- and state-based).


License
-------

Freely available subject to a BSD 2-Clause License.  
Please cite this code by adding the following reference to your work:

>[**[PMQ18] S. Pagani, A. Manzoni, A. Quarteroni. Numerical approximation of parametrized problems in cardiac electrophysiology by a local reduced basis method**, submitted to CMAME, 2018.]

Development
-------

LocalROM was developed and is currently maintained by [`Stefano Pagani`](https://stefanopagani.github.io).
