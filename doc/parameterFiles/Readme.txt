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

