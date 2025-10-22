% MLTED, ELTED, ZCTED, GTED, MMTED, GOTED , or GOTED_Mod  Lee ,  SLN
clc;close all;clear;
addpath("Phase_Sync\")
rolloff = 0.2;
plotSCurve = 1;
TED='Lee';
Kp_analytic_ml = calcTedKp(TED, rolloff, 'simulated', plotSCurve);