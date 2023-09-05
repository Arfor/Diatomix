classdef RbCs
    %KRB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        d0  = 1.225 * 3.33564e-30; %Dipole Moment (V/m), 1Debye = 3.33564e-30 C*m  from Ni et al., Science 322, 231-235 (2008) 0.574
    end
    properties
        name    %
        Atom1   %Rb
        Atom2   %Cs
        c1      %for Atom 1
        c2      %for Atom 2
        c3      %
        c4      %
        gr      %Molecular Magnetic g-factor
        Brot    %Rotational constant
        Drot    %Centrifugal constant
        Q1      %nuclear electric quadrupole constant (at position of K) 
        Q2      %nuclear electric quadrupole constant (at position of Rb) 
        a0      %h*Hz/(W/cm^2) at 1064nm
        a2      % tensor polarisability at 1064nm
    end
    
    methods
        function obj = RbCs(nRb,nCs)
            arguments
                nRb {mustBeMember(nRb,[87])} = 87;
                nCs {mustBeMember(nCs,[133])} = 133;
            end
            %RbCs Construct an instance of this class            
            if nCs==133&&nRb==87
                obj.Atom1=Rb87;
                obj.Atom2=Cs133;
                idx=1;
            end
            %   Constants taken from DOI: 10.1038/s41567-021-01328-7
            h = 6.62606896e-34; 
            c1 = [98.4] * h;%for Rb
            c2 = [194.2] * h;%for Cs
            c3 = [192.4] * h;
            c4 = [19.0189557e3] * h;
            gr = [0.0062];
            Brot = 1e9*[0.49174] * h; 
            Drot = [207.3] * h; 
            Q1 = 1e6*[-0.80929]*h; %nuclear electric quadrupole constant (at position of Atom 1) 
            Q2 = 1e3*[59.98]*h;%nuclear electric quadrupole constant (at position of Atom 2) 
            
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

