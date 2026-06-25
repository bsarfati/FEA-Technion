# FEA-Technion Repository

This repository contains all of my codes written for the Finite Element Analysis graduate course at the Technion (minimal to no AI used).

## Primary Contents

* `FEA_HW1.m`, `FEA_HW2.m` (actually HW3), `FEA_HW3.m` (actually HW4), and `FEA_HW5.m` are homeworks from the course.

* `Assigments/` contains the homework guidelines.

* `Submissions/` contains the homeworks' accompanying reports.
 
* `Old HWs/` contains assignments from a similar undergraduate course which I completed during my BSc.

* `finalMain.m`, `finalEigen.m`, and `finalErrorNorms.m` are the main scripts written for the final project, which included building a FE sovler.

* Many more helper scripts...

### Helper scripts

Many independently functional helper scripts were written to complete the course objectives:

* **Mesh Operations:** Scripts for element handling and mesh merging such as `innerMesh.m`, `meshElement.m`, `mergeMesh.m`, `mergeMesh_main.m`, and `mergeMeshSoln.m`.

* **Calculations & Quadrature:** Mathematical utility scripts like `gaussQuadratureTri.m`, `getLagrangePoly.m`, and `integrateGeom.m`.

* **Mapping & Plotting:** Visualizations and boundary mapping tools, including `mapLocalSolution.m`, `retrieveMapping.m`, `getOrderedBoundary.m`, and `plotElement.m`.

* In depth descriptions of each script are provided mainly in the **final project document** but also in **homework documents**.
