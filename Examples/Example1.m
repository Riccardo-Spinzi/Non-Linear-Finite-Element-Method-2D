function INPUT = Example1
% This example simulates a generic template to run the finite element code
% with the FEM_RUNME.m script. This particular problem is about a plate
% with a hole in the middle, pullled by its left and right sides. Only a
% quarter of the problem is studied to diminish the number of elements due
% to simmetry of the problem itself. The problem is studied as linear.
%
%  ------------- FOR MORE INFO READ THE "FEM_NL_Setup.txt" file

% -- Init
INPUT = struct();

% -- Elements
INPUT.integration_pts = 2;
INPUT.Barlow_points = 0;

INPUT.elements = [  1      6      7      2
 2      7      8      3
 3      8     10      4
 4     10     11      5
 6      9     12      7
 7     12     14      8
 9     13     16     12
 8     14     18     10
10     18     22     11
12     16     23     14
13     15     24     16
24     15     17     25
14     23     29     18
25     17     19     26
16     24     30     23
26     19     20     27
27     20     21     28
18     29     35     22
30     24     25     31
31     25     26     32
23     30     36     29
32     26     27     33
33     27     28     34
36     30     31     37
37     31     32     38
38     32     33     39
29     36     45     35
39     33     34     40
45     36     37     44
44     37     38     43
43     38     39     42
42     39     40     41
41     46     47     42
42     47     48     43
43     48     49     44
44     49     50     45
46     51     52     47
47     52     53     48
48     53     54     49
49     54     55     50
51     56     57     52
52     57     58     53
53     58     59     54
54     59     60     55 ];

% -- Nodes
INPUT.nodes = [        0.         25. 
       0.       31.25 
       0.        37.5 
       0.       43.75 
       0.         50. 
4.8772581   24.519632 
6.7829435   30.889724 
 8.688629   37.259816 
9.5670858   23.096988 
10.594315   43.629908 
     12.5         50. 
13.425314   29.822741 
13.889256    20.78674 
17.283543   36.548494 
 17.67767    17.67767 
19.791942   28.090055 
 20.78674   13.889256 
21.141771   43.274247 
23.096988   9.5670858 
24.519632   4.8772581 
      25.          0. 
      25.         50. 
25.694628    35.39337 
25.758252   25.758252 
28.090055   19.791942 
29.822741   13.425314 
30.889724   6.7829435 
    31.25          0. 
31.597314   42.696685 
33.838835   33.838835 
 35.39337   25.694628 
36.548494   17.283543 
37.259816    8.688629 
     37.5          0. 
     37.5         50. 
41.919417   41.919417 
42.696685   31.597314 
43.274247   21.141771 
43.629908   10.594315 
    43.75          0. 
      50.          0. 
      50.        12.5 
      50.         25. 
      50.        37.5 
      50.         50. 
58.333333          0. 
58.333333        12.5 
58.333333         25. 
58.333333        37.5 
58.333333         50. 
66.666667          0. 
66.666667        12.5 
66.666667         25. 
66.666667        37.5 
66.666667         50. 
      75.          0. 
      75.        12.5 
      75.         25. 
      75.        37.5 
      75.         50.  ];

% -- Section properties
INPUT.mat_type = 'svk'; 
INPUT.mat_prop  = [72000, .3,  1.];     % [E, nu, t]

% -- Loading conditions
INPUT.load = [   56  1 125.
 57  1 250.
 58  1 250.
 59  1 250.
 60  1 125. ];

% -- Boundary conditions
INPUT.spc = [  1    1
 2    1
 3    1
 4    1
 5    1
21    2
28    2
34    2
40    2
41    2
46    2
51    2
56    2 ];

% --- Nonlinear solver parameters 
INPUT.sol_type     = 'linear';
INPUT.formul       = 'UL';
INPUT.NR_param  = [1.0, 0.1, 1e-8, 100]; % lambda_max, dlambda, norm_res_max, nitermax
