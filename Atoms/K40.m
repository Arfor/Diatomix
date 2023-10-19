classdef K40 < Atom
%K40

properties (Constant)
  mass = 39.96399817*1.66053906660e-27; %kg
  name = "{}^{40}K";
  parity = 1; %fermion = 1, boson = 0
  spin = 4; 
  electricQuadrupoleMoment = 1e6*0.45 * 6.62606896e-34;%(eQq)_K (Hz) from Ospelkaus et al., PRL 104, 030402 (2010)
  gFactor = -0.324;%from 10.1103/PhysRevA.78.033434#
  nuclearShielding = 1321*1e-6;%sigma from 10.1103/PhysRevA.78.033434
end
properties
   transitions
end

methods
   function obj = K40()
        %            transition(1).Wavelength = 780.241*1e-9; %m
        % transition(1).Frequency = 384.2304844685e12; %Hz
        % transition(1).Gamma = 38.11*1e6;
        % transition(1).Isat = 560; %W/m^2
        % transitions = struct2table(transition);
                
        hbar = 1.0546e-34;
        kB = 1.3806e-23;
        c = 299792458;
        % transitions.Wavelength = c./transitions.Frequency;
        % transitions.WavevectorK= 2*pi./transitions.Wavelength; 
        % transitions.Isat = pi/3 * (2*pi*hbar*c*transitions.Gamma) ./ (transitions.Wavelength.^3);
        % transitions.TDoppler = hbar*transitions.Gamma/(2*kB);
        % transitions.TRecoil = hbar^2.*(transitions.WavevectorK.^2)/(obj.mass*kB);
        % obj.transitions = transitions;
   end 
end
end

