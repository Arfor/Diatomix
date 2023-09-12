classdef Rb87 
%RB87 Summary of this class goes here
%   Detailed explanation goes here

properties (Constant)
  mass = 86.909180527*1.66053906660e-27; %kg
  name = "{}^{87}Rb";
  parity = 0; %fermion = 1, boson = 0
  spin = 3/2;
  electricQuadrupoleMoment = 1e6*-1.41 * 6.62606896e-34;%(eQq)_K (Hz) from Ospelkaus et al., PRL 104, 030402 (2010)
  % gFactor = 1.834;%from 10.1103/PhysRevA.78.033434
  % nuclearShielding = 3469*1e-6;%sigma(ppm) from 10.1103/PhysRevA.78.033434

  gFactor = 1.8295;%for comparison with diatomic-py
  nuclearShielding = 0;
end
properties
   transitions
end

methods
   function obj = Rb87()
        %            transition(1).Wavelength = 780.241*1e-9; %m
        transition(1).Frequency = 384.2304844685e12; %Hz
        transition(1).Gamma = 38.11*1e6;
        transition(1).Isat = 560; %W/m^2
        transitions = struct2table(transition);
        
%         transitions(2,:) = array2table([377.1074635e12,1.0e6,1.1]);
        
        hbar = 1.0546e-34;
        kB = 1.3806e-23;
        c = 299792458;
        transitions.Wavelength = c./transitions.Frequency;
        transitions.WavevectorK= 2*pi./transitions.Wavelength; 
        transitions.Isat = pi/3 * (2*pi*hbar*c*transitions.Gamma) ./ (transitions.Wavelength.^3);
        transitions.TDoppler = hbar*transitions.Gamma/(2*kB);
        transitions.TRecoil = hbar^2.*(transitions.WavevectorK.^2)/(obj.mass*kB);
        obj.transitions = transitions;
   end 
end
end

