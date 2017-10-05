# How to use Fundamental Zone matlab code

In the FZ_routines folder there are two matlab functions: fz_inserter.m and mackenzie_plot.m.

Example usage of fz_inserter.m: 
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
