classdef Cs133 
%RB87 Summary of this class goes here
%   Parameters from Aldegunde & Hutson (2017): 10.1103/PhysRevA.96.042506

properties (Constant)
  mass = 132.905451933*1.66053906660e-27; %kg
  name = "{}^{133}Cs";
  parity = 0; %fermion = 1, boson = 0
  spin = 7/2;
  gFactor = 0.738;%
  % electricQuadrupoleMoment = 1e6*-1.41 * 6.62606896e-34;%(eQq)_Cs (Hz) 
  nuclearShielding = 6278.7e-6;%sigma(ppm)

  % gFactor = 0.7331; %for comparison with diatomic-py
  % nuclearShielding = 0;
end
properties
   transitions
end

methods
   function obj = Cs133()
        %            transition(1).Wavelength = 780.241*1e-9; %m
%         transition(1).Frequency = 384.2304844685e12; %Hz
%         transition(1).Gamma = 38.11*1e6;
%         transition(1).Isat = 560; %W/m^2
%         transitions = struct2table(transition);
% 
% %         transitions(2,:) = array2table([377.1074635e12,1.0e6,1.1]);
% 
%         hbar = 1.0546e-34;
%         kB = 1.3806e-23;
%         c = 299792458;
%         transitions.Wavelength = c./transitions.Frequency;
%         transitions.WavevectorK= 2*pi./transitions.Wavelength; 
%         transitions.Isat = pi/3 * (2*pi*hbar*c*transitions.Gamma) ./ (transitions.Wavelength.^3);
%         transitions.TDoppler = hbar*transitions.Gamma/(2*kB);
%         transitions.TRecoil = hbar^2.*(transitions.WavevectorK.^2)/(obj.mass*kB);
%         obj.transitions = transitions;
   end 
end
end

