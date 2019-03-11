#Script for PS3-S19

#Magic function
include("Flux.jl");

#Plot plot plot
using PyPlot

#Parameters
E = 0.01/1000; #steady-state enzyme concentration [mmol/gDW]
b1 = 10/3600; #Carbanoyl Phosphate input [mmol/gDW s]
b2 = (10/3600)*((1.49*10^(-2))/(1.49*10^(-2)+1*10^(-2))); #Aspartate input [mmol/gDW s]
b3 = (10/3600)*((4.85*10^(-4))/(4.85*10^(-4)+5.3*10^(-3))); #Fumarate output [mmol/gDW s]
b4 = 10/3600; #Urea output [mmol/gDW s]
b5 = (10/3600)*((4.67*10^(-3))/(4.67*10^(-3)+3*10^(-5))); #ATP input [mmol/gDW s]
b6 = 10/3600; #Ammonia input [mmol/gDW s]
b7 = 10/3600; #Orthophosphate output [mmol/gDW s]
b8 = (10/3600)*((4.23*10^(-5))/(4.23*10^(-5)+6.46*10^(-5))); #AMP output [mmol/gDW s]
b9 = 10/3600; #Diphosphate output [mmol/gDW s]
v1 = 203*E; #6.3.4.5 [mmol/gDW s]
v2 = 34.5*E; #4.3.2.1 [mmol/gDW s]
v3 = 249*E; #3.5.3.1 [mmol/gDW s]
v4 = 88.1*E; #2.1.3.3 [mmol/gDW s]
v5 = (13.7*E)*((2.55*10^(-4))/(2.55*10^(-4)+2.5*10^(-3))); #1.14.13.39 [mmol/gDW s]

td = 20*3600; #Doubling time [mmol/gDW s]

#Stochiometric matrix.
#Columns: v1,v2,v3,v4,v5,b1,b2,b3,b4,b5,b6,b7,b8,b9
#rows: C,Ort,O,CP,ATP,A,AMP,PPi,Arsu,F,Arg,U,H2O,NH3
stoichiometric_matrix = [[-1.0 0 0 1 -1 0 0 0 0 0 0 0 0 0];
                         [0 0 0 1 0 0 0 0 0 0 0 -1 0 0];
                         [0 0 1 -1 0 0 0 0 0 0 0 0 0 0];
                         [0 0 0 -1 0 1 0 0 0 0 0 0 0 0];
                         [-1 0 0 0 0 0 0 0 0 1 0 0 0 0];
                         [-1 0 0 0 0 0 1 0 0 0 0 0 0 0];
                         [1 0 0 0 0 0 0 0 0 0 0 0 -1 0];
                         [1 0 0 0 0 0 0 0 0 0 0 0 0 -1];
                         [1 -1 0 0 0 0 0 0 0 0 0 0 0 0];
                         [0 1 0 0 0 0 0 -1 0 0 0 0 0 0];
                         [0 1 -1 0 1 0 0 0 0 0 0 0 0 0];
                         [0 0 1 0 0 0 0 0 -1 0 0 0 0 0];
                         [0 0 -1 0 1 0 0 0 0 0 0 0 0 0];
                         [0 0 0 0 -1 0 0 0 0 0 1 0 0 0]];


Lb = zeros(Int8, 8, 1); #Lower bound array
Ub = [10;10;10;10;10;10;10;10]; #Upper bound array

#Bounds array
default_bounds_array = [[0 v1];
                        [0 v2];
                        [0 v3];
                        [0 v4];
                        [-v5 v5];
                        [0 b1];
                        [0 b2];
                        [0 b3];
                        [0 b4];
                        [0 b5];
                        [0 b6];
                        [0 b7];
                        [0 b8];
                        [0 b9]];
#Species bounds array
species_bounds_array = zeros(Float64, 14, 2);
#Objective array
objective_coefficient_array = [0.0; 0; 0; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0]; #Change -1 to max or 1 to min

#Flux calculation
flux, n1, n2, n3, n4, n5 = calculate_optimal_flux_distribution(stoichiometric_matrix,default_bounds_array,species_bounds_array,objective_coefficient_array);
flux*3600*1000 #Flux in [umol/gDW h]
