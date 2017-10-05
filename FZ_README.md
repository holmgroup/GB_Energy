# How to use Fundamental Zone matlab code

In the FZ_routines folder there are two matlab functions: fz_inserter.m and mackenzie_plot.m.

```Matlab
function [aa_pairs_fz, disorientation_fz, OM_fz, r_fz, th_fz, sizes] = fz_inserter(O_list1,O_list2,plotting_str)
```
The input for the code is two N x 9 lists of orientations O_list1, O_list2. 
For the sake of example, let N = 1: 

Define two orientation matrices O1 and O2, with crystallographic directions [h k l] as rows. 

```Matlab
O1 = [1 1 1; 1 -1 0; 1 1 -2];
O2 = [1 1 0; 1 -1 0; 0 0 1];

O1f = reshape(O1',[1 9]);
O2f = reshape(O2',[1 9]);

[aa_example,~,~,~,~,~] = fz_inserter(O1f,O2f);
```
fz_inserter takes flattened orientation matrices in 1 x 9 format [ row1 row2 row3 ], represented by O1f and O2f above, and finds the rotation between those orientations with minimum rotation angle over all cubic symmetry operators (for angular component of axis angle pair). The rotation axis of the axis angle pair is placed in the standard stereographic triangle (SST), of which there are 24 arbitrary choices. This code calculates the fundamental rotations between rows of two Nx9 input matrices. See matlab function for more information on output options. 

2.Example usage of mackenzie_plot.m:
```Matlab
function [disorientation_list, axis_list] = mackenzie_plot(n,input_str)
%example usage, generate mackenzie distribution with 1000 random rotations 
[dlist, alist] = mackenzie_plot(1000,'on'); 
%dlist is list of disorientation angles
%alist is list of axes placed in SST
```
