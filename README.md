<p align=center><img height="100%" width="100%" src="images/logo.png"></p>

_Anura3D_ is a software for the numerical modelling and simulation of large deformations and
soil–water–structure interaction using the material point method (MPM). Copyright (C) 2020 Members of the Anura3D MPM Research Community.

**Anura3D** is **free** you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation <https://www.gnu.org/licenses/>, either version 3 of the License, or (at your option) any later version

Anura3D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.


# Main Features

The features implemented in the open source are briefly summarized in the Table below. While some of them have been tested, others are still under development. The examples provided in the tutorial provide guidance to the user through the available tested features.

_State legend_:

✔️: Available and tested

🚧: Under development

|Feature| State|
|---|---|
|**Geometrical dimensions**| |
|2D-plane strain | ✔️|
|2D-axisymmetric | 🚧|
|3D-cartesian| ✔️|
|3D-cylindrical| ✔️|

The current state of the Anura3D has the following limitations concerning the tested features:


* 3D-Cylindrical: y-axis is the axis of symmetry, gravity can only be applied in y-axis direction.
* Multiphase formulations: a combination of different material types can cause problems for certain combinations (undrained total stress and undrained effective stress), no water flow is transferred between saturated and dry materials.
* Contact algorithm: the maximum number of contact materials is four (4), the maximum number of master materials is one (1). Contact is fixed on the mesh nodes.
* Moving mesh: works only for prismatic bodies in 3D and trapezium areas in 2D, the moving mesh direction is constrained only in one direction. the moving mesh needs an extension and/or compression mesh.
* Excavation: limited to 30 excavation stages.
* Rigid body: can be applied only to one body in the system, only works together with the contact algorithm, only moves in one direction, it can’t rotate.
* Fixities and traction boundary conditions: only aligned with the coordinate axes
* K0 stress initialization: limited only to horizontal surfaces, one homogeneous material (otherwise, use gravity together with local damping and stress initialization with quasistatic convergence criteria).
The user should also take into account the following warnings:
* Multiple materials: the use of several drainage material types on the same model is not recommended.
* Absorbing boundaries: do not prevent material points from leaving the mesh.



# Documentation
Here you can find how to get started with _Anura3D_:
## Tutorial manual
A comprehensive guide that shows how to compile, preprocess, and postprocess problems with Anura3D, GID, and Paraview.

* Link to the manual.

You can also watch our tutorials in our YouTube channel. A link to the different tutorial videos is below.

* [Introduction to Anura3D](https://www.youtube.com/watch?v=6Rx98oyO51A)
* [How to download and compile Anura3D](https://www.youtube.com/watch?v=1qlRcZvAZ_A)
* [Overview of the calculation process](https://www.youtube.com/watch?v=-kbWmlQrfao)
* [One-dimensional consolidation (oedometric compression test)](https://www.youtube.com/watch?v=nvOIf05ie4k)
* [Triaxial compression test](https://www.youtube.com/watch?v=a7qNT_Qo8Tg)
* [Sliding blocks] 🚧
* [Column collapse](https://www.youtube.com/watch?v=3dvcvYgI2cIs)
* [Shallow Foundation](https://www.youtube.com/watch?v=p4YjcpIk0uE)
* Impact problem [2D](https://www.youtube.com/watch?v=qiXNbpx0vZc) and [3D](https://www.youtube.com/watch?v=fNlEdc_nhgo)
* [Excavation](https://www.youtube.com/watch?v=TbTj6upumqs)
* [Submerged slope collapse](https://www.youtube.com/watch?v=kucDI3AKnRY)
* [Shallow foundation in 3D-cylindrical coordinates] 🚧


## Scientific manual

Learn about the theoretical models and numerical techniques implemented in _Anura3D_.

* Link to the manual

## Coding manual

We love people that modifies our code and develops useful and cool features. We provide the _coding manual_ to ensure uniformity of coding style, documentability, generality, and compatibility with existing features.

* link to the manual


# Examples of use
_Anura3D_ has been extensively used in various problems in geotechnical and coastal engineering. Below you can find videos of the usage of _Anura3D_ and links to presentations of members of the Anura3D Research Community.

[![Coseismic landslides and internal erosion](https://img.youtube.com/vi/Cd54tmVGG84/0.jpg)](https://www.youtube.com/watch?v=Cd54tmVGG84)
[![An unsaturated formulation in MPM for dams and levees](https://img.youtube.com/vi/K0zl_Q3S6uM/0.jpg)](https://www.youtube.com/watch?v=K0zl_Q3S6uM)
[![MPM simulation of penetration problems](https://img.youtube.com/vi/RHdNfkqyYqQ/0.jpg)](https://www.youtube.com/watch?v=RHdNfkqyYqQ)

# News and events

Stay tuned.

# Anura3D MPM Research Community

# Publications

# Special Thanks
The support of CIMNE <https://www.cimne.com/> is greatly acknowledged.

# Users
Anura3D is used by:

# Get in touch

[info@anura3D.com](mailto:info@anura3d.com)
