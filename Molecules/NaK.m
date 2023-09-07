classdef NaK
    %KRB Summary of this class goes here
% Parameters from Will et al., PRL 116, 225306 (2016)
% https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.225306
% and from Aldegunde et al., PRA 96, 042506 (2017)
% https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.042506
    
    properties (Constant)
        d0  = 2.72 * 3.33564e-30; %Dipole Moment (V/m), 1Debye = 3.33564e-30 C*m 
    end
    properties
        name    %
        Atom1   %Na
        Atom2   %K
        c1      %for Atom 1
        c2      %for Atom 2
        c3      %
        c4      %
        gr      %Molecular Magnetic g-factor
        Brot    %Rotational constant
        Drot    %Centrifugal constant
        Q1      %nuclear electric quadrupole constant (at position of Atom1) 
        Q2      %nuclear electric quadrupole constant (at position of Atom2) 
        a0      %h*Hz/(W/cm^2)
        a2      % tensor polarisability
    end
    
    methods
        function obj = NaK(nNa,nK)
            arguments
                nNa {mustBeMember(nNa,[23])} = 23;
                nK {mustBeMember(nK,[40])} = 40;
            end
            %RbCs Construct an instance of this class            
            if nK==40&&nNa==23
                obj.Atom1=Na23;
                obj.Atom2=K40;
                idx=1;
            end
            %   Constants taken from DOI: 10.1038/s41567-021-01328-7
            h = 6.62606896e-34; 
            c1 = [117.4] * h;%for Atom1
            c2 = [-97.0] * h;%for Atom2
            c3 = [-48.4] * h;
            c4 = [-409] * h;
            gr = [0.0253];
            Brot = 1e9*[2.8217297] * h; 
            Drot = [207.3] * h; 
            Q1 = 1e6*[-0.187]*h; %nuclear electric quadrupole constant (at position of Atom 1) 
            Q2 = 1e6*[0.899]*h;%nuclear electric quadrupole constant (at position of Atom 2) 
            
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

