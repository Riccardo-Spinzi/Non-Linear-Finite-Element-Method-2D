%% F.E. METHOD
clear
close all
clc

addpath('C:\Users\Utente\Desktop\FEM NL\INPUT files');

% ---1. Pre-process

INPUT = inputNL_cantilever_5_els('neohookean','ul');          

% if the problem is non-linear, specify the Newton-Rhapson method to apply:
NR_method = 'NR';

% ---2. Solution

[ MODEL, SOL, POST ] = fem_solver( INPUT, NR_method );

% ---3. Post-process ( no 3D ) 

plot_shapes( MODEL );

[ MODEL ] = stress_plot( MODEL, POST, SOL );

plot_eq_path( MODEL, POST, SOL )

Display_results( MODEL, SOL, POST )















