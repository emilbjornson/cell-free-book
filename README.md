[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=emilbjornson/cell-free-book)

Foundations of User-Centric Cell-Free Massive MIMO
==================

This repository contains the free authors' version and the accompanying code package of the textbook :

Özlem Tuğfe Demir, Emil Björnson, and Luca Sanguinetti (2021) “[Foundations of User-Centric Cell-Free Massive MIMO](https://www.nowpublishers.com/article/Details/SIG-109)”, Foundations and Trends in Signal Processing: Vol. 14, No. 3-4, pp. 162-472. DOI: 10.1561/2000000109.

The code package contains a simulation environment, based on Matlab, that can be used to reproduce all the simulation results in the monograph. We hope that the code will support you in the learning of the Cell-free Massive MIMO topic and also serve as a baseline for further research endeavors. *We encourage you to also perform reproducible research!*

## Abstract of the Book

Imagine a coverage area where each mobile device is communicating with a preferred set of wireless access points
(among many) that are selected based on its needs and cooperate to jointly serve it, instead of creating autonomous
cells. This effectively leads to a user-centric post-cellular network architecture, which can resolve many of the interference issues and service-quality variations that appear in cellular networks. This concept is called User-centric Cellfree Massive MIMO (multiple-input multiple-output) and
has its roots in the intersection between three technology
components: Massive MIMO, coordinated multipoint processing, and ultra-dense networks. The main challenge is to
achieve the benefits of cell-free operation in a practically
feasible way, with computational complexity and fronthaul
requirements that are scalable to enable massively large
networks with many mobile devices. This monograph covers
the foundations of User-centric Cell-free Massive MIMO,
starting from the motivation and mathematical definition. It
continues by describing the state-of-the-art signal processing
algorithms for channel estimation, uplink data reception, and downlink data transmission with either centralized or
distributed implementation. The achievable spectral efficiency is mathematically derived and evaluated numerically
using a running example that exposes the impact of various
system parameters and algorithmic choices. The fundamental tradeoffs between communication performance, computational complexity, and fronthaul signaling requirements
are thoroughly analyzed. Finally, the basic algorithms for
pilot assignment, dynamic cooperation cluster formation,
and power optimization are provided, while open problems
related to these and other resource allocation problems are
reviewed. All the numerical examples can be reproduced
using the accompanying Matlab code.

## Content of Code Package

This code package contains 30 Matlab scripts and 17 Matlab functions. You can run the code in MATLAB online without a license by clicking on the link above.

Each script is used to reproduce a particular simulation-generated figure in the book. The scripts are named using the convention chapterX_figureY, which is interpreted as the script that reproduces Figure X.Y. A few scripts are instead named as sectionX_figureY_Z and will then generate both Figure X.Y and Figure X.Z.

The functions are used by the scripts to carry out certain tasks, such as initiating a simulation setup, generating channel correlation matrices, generating channel realizations, computing channel estimates, computing SEs, implementing power control algorithms, etc.

See each script and function for further documentation. Note that some of the functions use [CVX](http://cvxr.com/cvx/) with SDPT3 solver, which need to be installed separately.

## Acknowledgements

We would first like to thank our students and collaborators in the areas
related to this monograph. Without the results, encouragements, and
insights obtained through our joint research during the last decade, it
wouldn’t have been possible to write this monograph. We are grateful for
the constructive feedback from the reviewers, which helped us to focus
our final editing efforts at the right places. In particular, we would like
to thank Angel Lozano, Jiayi Zhang, Mahmoud Zaher, and Yasaman
Khorsandmanesh for giving detailed comments.

Özlem Tuğfe Demir and Emil Björnson have been supported by the
Wallenberg AI, Autonomous Systems and Software Program (WASP)
funded by the Knut and Alice Wallenberg Foundation. Emil Björnson
has also been supported by the Excellence Center at Linköping – Lund in
Information Technology (ELLIIT), the Center for Industrial Information
Technology (CENIIT), the Swedish Research Council, and the Swedish
Foundation for Strategic Research. Luca Sanguinetti has been partially
supported by the University of Pisa under the PRA Research Project
CONCEPT, and by the Italian Ministry of Education and Research
(MIUR) in the framework of the CrossLab project (Departments of
Excellence).

## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our textbook as described above. We also recommend that you mention the existence of this code package in your manuscript, to spread the word about its existence and to ensure that you will not be accused of plagiarism by the reviewers of your manuscript.
