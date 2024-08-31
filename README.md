## ZStruct
code for ZStruct-2

## Overview

ZStruct is a package for automated reaction discovery developed in C++.

For more information, check out the wiki page (somewhat out of date but hopefully still useful): https://github.com/ZimmermanGroup/ZStruct/wiki

Sample tutorial files can be found under the tutorial folder: https://github.com/ZimmermanGroup/ZStruct/tree/master/tutorial

## Installation

First, clone this github repository.

Aug 2024: Attempting a more reproducible build with [Pixi](https://pixi.sh), a tool built on the conda ecosystem to manage both dependencies and configuration / installation workflows.

The pixi.toml file in this repository specifies the dependency package names as they appear on [conda-forge](https://conda-forge.org/) and installation tasks. You can use pixi (simple installation instructions [here](https://pixi.sh)) to install the dependencies and execute the build, or reference the configuration files and do something equivalent using your preferred tools.

With pixi, I built ZStruct by navigating into the ZStruct git repository and executing `pixi run start` which installs the dependencies, executes the Makefile, and if successful runs the built executable.

Notes: The build has succeeded in a github actions environment which is a freshly cloned linux virtual machine. On our computing cluster, I find that I either need to unload one of our default linux modules or clear the LD_LIBRARY_PATH with `export LD_LIBRARY_PATH=""`. I assume this avoids a conflicting version of some dependency taking priority over the version installed for ZStruct.
