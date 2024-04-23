# Non-Linear-Finite-Element-Method
F.E.M. repository ~ WIP

This code comes as a really basic algorithm to evaluate displacements and forces in a 2D structure 
modeled with linear and non linear continuum elements. To interpolate the results Gauss integration is used, 
defining integration points inside a single element as checkpoints to evaluate the results. 
To use the program open the "my_struct_name".m file and write the specific structure's properties, following
the Example1.m file. To have a complete overview of the program, please read the FEM_NL_Setup.txt file.
To modify the scale of the displacements figures open the plot_shapes.m file and change the 
n_visual parameter. 

Implemented features:

~ Solve problems with 2D non linear continuum elements;

~ Get strains (only NL) and stresses in the structure;

~ Plot the equilibrium path to visualize the non linear behavior of the structure;

~ Plot undeformed and deformed structures at the same time to compare them; 

~ Change between pure Newton-Rhapson, modified Newton-Rhapson and initial stiffness evaluation methods;

~ Evaluate different kinds of hyperelastic constitutive laws.

Future aims:

~ Implement the possibility to evaluate distributed forces along the structure;

~ Implement strain recovery for the linear case;

~ Implement prescribed displacements as boundary conditions in the model;
 
~ Implement a 3D version of the algorithm, for both linear and non linear continuum elements;

~ Implement beam elements for a generic non linear structure;

~ Implement simple plasticity models to account for linear plasticity;

~ Implement a displacement based approach to evaluate the non linear solution;

~ Implement path following techniques to solve highly non linear problems and avoiding NR singularity;

~ Implement a python code to substitute the now published one and extend it to non-matlab users. 


NOTE: future aims are not to be intended in the order they are presented. Depending on what I think I could
implement first, it will be done.

Changelog:

----23-04 RELEASED v1.0

First draft of the program released, along with full documentation, flowchart and examples.
