**** Description of parameter Files used inside FEDDLib ****
@2024

There are several parameter files used in the FEDDLib that play a crucial role in defining and solving physical problems. These parameter files allow the user to customize various aspects of the problem and solution parameters.
Here is a breakdown of the three main parameter files:

1. parametersProblem.xml: This file contains all the parameters related to the specific details of the definition of a problem. It includes physical constants, like e.g. density, or solver related information, like e.g. solver tolerances, as well as options for visualization export. 
                          By modifying this file, the user can adjust the problem-specific settings.

2. parametersSolver.xml:  This file includes parameters that can be used to specify the setting of the linear solver and nonlinear solver. Here, the linear solver type is often Belos, which indicates the utilization of the Belos package. 
                          If the user sets in the parametersProblem.xml file, the linearization to be NOX they can specify here NOX-related parameters, such as globalization strategies and related parameters.

3. parametersPrec.xml:    This file is dedicated to parameters related to the usage of preconditioners. If the parametersProblem.xml file specifies the Monolithic option, the preconditioner type should be set to Frosch, as it is implemented in the Frosch package. 
                          On the other hand, if the Blockpreconditioner Teko was specified, the corresponding Teko-related parameters should be set in this file.

We include here four different files that provide an overview of the meaning and options of different parameters inside these files. The  parametersPrec.xml File is here further divided into parametersPrecMonolithic.xml and parametersPrecTeko.xml files, as they have distinct keywords. 
It's important to note that this overview is not intended to be copied directly but rather to provide a comprehensive understanding of the parameter files and their purpose. If any new relevant keywords are introduced, they should be added to the list along with an explanation of their usage and significance.


-----------------------------------------------------------------------------

Collection of detailed explanations
===================================


-----------------------------------
RGDSW coarse space -- Coarse space options

RGDSW: reduced-dimension GDSW
GDSW: generalized Dryja-Smith-Widlund

The implemented coarse spaces of the RGDSWCoarseOperator are based on the following paper:

(*) Clark R. Dohrmann, Olof B. Widlund, 2017, On the Design of Small Coarse Spaces for Domain Decomposition Algorithms
    https://doi.org/10.1137/17M1114272

Dohrmann/Widlund have proposed several variants of energy-minimizing coarse spaces. 
(Technically, they are only energy-minimizing for symmetric, positive definite problems.) 
These prescribe a partition of unity on the interface, which is subsequently extended to the interior of the subdomains.
Not all variants are implemented in FROSch (which FEDDLib uses). 
Furthermore, the variant coming from adaptive coarse spaces is also not implemented; cf. Fig. 1 in
    A. Heinlein, A. Klawonn, J. Knepper, O. Rheinbach, O. B. Widlund, 2022, Adaptive GDSW Coarse Spaces of Reduced Dimension for Overlapping Schwarz Methods
    https://doi.org/10.1137/20M1364540

The variants of the RGDSW coarse spaces in (*) are named Option 1, Option 2.1, Option 2.2, and Option 2.1+2. 
Thereof, Option 1 and Option 2.2 are implemented.

Option 1: Multiplicity scaling on the interface
   Prescribe the inverse of the "coarse node multiplicity" to a node/dof on the interface.
   Eq. (1) in (*)
   Coarse node multiplicity: Let n be a node/dof. How many coarse nodes are ancestors of n?
      Ancestor: The set of adjacent subdomains of a node/dof is a subset of (not equal to) the set of adjacent subdomains of the ancestor.
      Detailed description: RGDSW generates interface components that overlap.
      Each interface component is (generally) associated with a coarse node. 
      For each coarse node, adjacent coarse edges and coarse faces are added to obtain a "coarse star".
      This coarse star is the RGDSW interface component and overlaps with neighboring coarse stars.
      The question now is: To which coarse stars does a node belong? This is the sought multiplicity.

Option 2.2: Distance scaling on the interface
   Prescribe sth. like (1/d1(n)) / (1/d1(n) + ... + 1/dk(n)) to a node/dof n.
   di(n) is the euclidian distance of n to its ancestor coarse node i in {1,...,k}.
   Eq. (5) in (*)


Example for choosing the coarse space variants in the parameter file:

<Parameter name="CoarseOperator Type" type="string" value="RGDSWCoarseOperator"/>
<ParameterList name="RGDSWCoarseOperator">
    <ParameterList name="Blocks">
        <ParameterList name="1">
            <Parameter name="Option" type="string" value="1"/>
            or 
            <Parameter name="Option" type="string" value="2.2"/>

-----------------------------------
