# WARR
Inferring density, radio-wave velocity and reflector depths from radar wide-angle surveys  

% Description:  
%----------------------------------------------------------------------------------------------  
% This code exemplifies how to invert radar-wide angle data using raytracing and Gauss-Newton.   
% The shape of the density profile is prescribed as a function of:  
%  
% density(depth) = 910 - A*exp(-r*depth)  
%  
% and the parameters A,r are inverted for together with reflector depths. This is a synthetic examples,    
% it can quickly be adapted for real data (real-case example should appear soon, contact me if not).  
% Refer to the following paper for details (citations are appreciated):  
%  
% "Density anomaly in ice-shelf channels as seen by optical televiewing and wide-angle radar"  
%  R. Drews, J. Brown, K. Matsuoka, E. Witrant, M. Philippe, B. Hubbard, F. Pattyn;   
%  submitted to The Cryosphere Discussions  
%  
%  
% Code has been tested on Matlab 2011b.  
%----------------------------------------------------------------------------------------------
  
  
% Requirements:     
% ---------------------------------------------------------------------------------------------  
% The raytracing forward model is based on     
% Margrave, G. F.: Numerical Methods of Exploration Seismology with algorithms in Matlab, CREWES    
% Toolbox Version: 1006 and requires the following files in a local folder "raytracing":     
%  
% drawray.m  flipy.m  shootray.m  sphdiv.m  surround.m  traceray_pp.m  
%  
% These can be downloaded here: http://www.codeforge.com/article/174599  
%----------------------------------------------------------------------------------------------  
   
  
% To Do:   
%---------------------------------------------------------------------------------------------  
% 1. Get Graphical User Interface to make it more accessible  
% 2. Use ASCII txt files for input data and unify scripts for real and syntethic cases  
% 3. Ramping of Layers (to simulate systematic bias, has currently hardcoded shape of slope.)  
% 3. Work on better documentation.  
% 4. Transfer everything to Python (GUI as well?)   

% Last Changes: 28.09.2015   
%----------------------------------------------------------------------------------------------  

% How to Use: The script is documented in file and hopefully somewhat self-explanatory.  
