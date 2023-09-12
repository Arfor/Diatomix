classdef KRb
    %KRB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        % d0  = 0.566 * 3.33564e-30; %Dipole Moment (V/m), 1Debye = 3.33564e-30 C*m  from Ni et al., Science 322, 231-235 (2008) 0.574
        d0  = 0.573999 * 3.33564e-30; %Dipole Moment (V/m) Till
    end
    properties
        name    %
        Atom1   %K
        Atom2   %Rb
        c1      %for Atom 1
        c2      %for Atom 2
        c3      %
        c4      %
        gr      %Molecular Magnetic g-factor
        Brot    %Rotational constant
        Drot    %Centrifugal constant
        Q1      %nuclear electric quadrupole constant (at position of Atom1) 
        Q2      %nuclear electric quadrupole constant (at position of Atom2) 
        a0      %h*Hz/(W/cm^2) at 1064nm
        a2      % tensor polarisability at 1064nm
    end
    
    methods
        function obj = KRb(nK,nRb)
            arguments
                nK {mustBeMember(nK,[39,40,41])} = 40;
                nRb {mustBeMember(nRb,[85,87])} = 87;
            end
            %KRB Construct an instance of this class            
            if nK==39&&nRb==85
                obj.Atom1=K39;
                obj.Atom2=Rb85;
                idx=1;
            elseif nK==39&&nRb==87
                obj.Atom1=K39;
                obj.Atom2=Rb87;
                idx=2;
            elseif nK==40&&nRb==85
                obj.Atom1=K40;
                obj.Atom2=Rb85;
                idx=3;
            elseif nK==40&&nRb==87
                obj.Atom1=K40;
                obj.Atom2=Rb87;
                idx=4;
            elseif nK==41&&nRb==85
                obj.Atom1=K41;
                obj.Atom2=Rb85;
                idx=5;
            elseif nK==41&&nRb==87
                obj.Atom1=K41;
                obj.Atom2=Rb87;
                idx=6;
            end
            %   Constants taken from DOI: 10.1103/PhysRevA.78.033434
            % Using Till's constants for K40Rb87
            h = 6.62606896e-34; 
            c1 = [19.9,19.8,-24.2,-24.1,10.5,10.4] * h;%for K
            c2 = [127.0,427.5,124.8,420.1,122.8,413.1] * h;%for Rb
            c3 = [11.5, 38.9,-14.2,-48.2,6.3,21.3] * h;
            c4 = [482.5,1635.7,-599,-2030.4,264.3,896.2] * h;
            gr = [0.0144,0.0142,0.0141,0.0140,0.0139,0.0138];
            Brot = 1e9*[1.142,1.134,1.123,1.1139514,1.104,1.096] * h; 
            Drot = [0,0,0,0,0,0] * h; 
            % QK =    1e6*[-0.245,-0.245,0.306,0.306,-0.298,-0.298]*h; %nuclear electric quadrupole constant (at position of K) 
            % QRb =   1e6*[-3.142,-1.520,-3.142,-1.520,-3.142,-1.520]*h;%nuclear electric quadrupole constant (at position of Rb) 
            QK =    1e6*[-0.245,-0.245,0.452,0.452,-0.298,-0.298]*h; %Diatomic-py / Till's constants
            QRb =   1e6*[-3.142,-1.308,-3.142,-1.308,-3.142,-1.308]*h;
            a0 =   1e6*[0,0,0,0.553,0,0]*h; %h*Hz/(W/m^2) 
            a2 =   1e6*[0,0,0,0.447,0,0]*h;
            obj.c3 = c3(idx);
            obj.c4 = c4(idx);
            obj.gr = gr(idx);
            obj.Brot = Brot(idx);
            obj.Drot = Drot(idx);
            obj.c1 = c1(idx);
            obj.c2 = c2(idx);
            obj.Q1 = QK(idx);
            obj.Q2 = QRb(idx);
            obj.a0 = QRb(idx);
            obj.a2 = QRb(idx);

            obj.name = obj.Atom1.name + obj.Atom2.name;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

