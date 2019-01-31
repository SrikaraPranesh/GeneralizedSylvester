% This script generates the results in the paper 
% "Perturbation analysis of multi-term generalised 
% Sylvester equation and its application to stochastic 
% Galerkin method -- S. Pranesh"



clear all; close all

addpath('AnalysisTestScripts')
addpath('MainScripts')
addpath('StochasticGalerkinScripts')


% GSEResultsGenerate;
Sylvester_StochasticGalerkin1;
Sylvester_StochasticGalerkin2;


movefile('*.txt','results')


