%% Read me

% this folder contains all necessary script and functions for the Burgers
% euqation 1D, based on finite element, model order reduction (via POD)

% -There is indeed no physical problem under this equation.
% -Most of the time this equation is used in order to test a numerical schematics.
% -However this equation is well-kown to show shock formation.
% -It is the 1D Navier-Stokes equation without the pressure term and the volume forces.

%  Modification:
%  WeiX, May 28th 2015, First Edition.
%  WeiX, Spe 30th 2016, Second version.


% functionName_Demo is the demonstration file for function:functionName

% Phi  = FEM_HatFunc(i,N,x )
% The hat function. 0<=x<=1!!

% Phi  = FEM_HatFunc_Diff(i,N,x )
% First order differential of the hat function. 0<=x<=1!!

