classdef NaCs
    %NaCs Summary of this class goes here
%   Parameters from Aldegunde & Hutson (2017): 10.1103/PhysRevA.96.042506
    properties (Constant)
        d0  = 4.75 * 3.33564e-30; %Dipole Moment (V/m), 1Debye = 3.33564e-30 C*m 
    end
    properties
        name    %
        Atom1   %
        Atom2   %
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
        function obj = NaCs(nNa,nCs)
            arguments
                nNa {mustBeMember(nNa,[23])} = 23;
                nCs {mustBeMember(nCs,[133])} = 133;
            end
            %RbCs Construct an instance of this class            
            if nCs==133&&nNa==23
                obj.Atom1=Na23;
                obj.Atom2=Cs133;
                idx=1;
            end
            %   Constants taken from DOI: 
            h = 6.62606896e-34; 
            c1 = [14.2] * h;%for Na
            c2 = [854.5] * h;%for CS
            c3 = [105.6] * h;
            c4 = [3941.8] * h;
            gr = [0.0062];
            Brot = 1e9*[1.7396] * h; 
            Drot = [0] * h; 
            Q1 = 1e6*[-0.097]*h; %nuclear electric quadrupole constant (at position of Atom 1) 
            Q2 = 1e6*[0.15]*h;%nuclear electric quadrupole constant (at position of Atom 2) 

            obj.c1 = c1(idx);
            obj.c2 = c2(idx);            
            obj.c3 = c3(idx);
            obj.c4 = c4(idx);
            obj.gr = gr(idx);
            obj.Brot = Brot(idx);
            obj.Drot = Drot(idx);
            obj.Q1 = Q1(idx);
            obj.Q2 = Q2(idx);

            obj.name = obj.Atom1.name + obj.Atom2.name;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

