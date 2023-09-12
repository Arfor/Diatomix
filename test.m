clear all
addpath("Molecules\");addpath("Atoms\")
C = Constants;
%% Define Fields
maxN = 2;
B.value = 0; %549.%
B.dir = [0,0,1];
E.value = 0;
E.dir = [0,0,1];
I.value = 0;
I.dir = [0,0,1];
I.pol = [1,0,0];
Fields.B = B;
Fields.E = E;
Fields.I = I;
clear E I B
%% Define Hamiltonian
Mol = RbCs;
Ham = Hamiltonian(Molecule = Mol, maxN=maxN);

H0 = Ham.hyperfine.total;
% + B.value*Ham.zeeman + E.value*Ham.dc_stark;
% H_ac = Ham.ac_stark;
Base = Ham.Basis;
nStates = Base.NStates;

Bs = linspace(0,50,50)*1e-4;
zeemanMap = nan(length(Bs),nStates);
for k = 1:length(Bs)
    B = Bs(k);
    H = H0 + B*Ham.zeeman; %why this sqrt(2) necessary to coincide with diatomic-py?
    [V,D] = eig(full(H/C.h),'vector'); %no need to use sparse matrices for matrices smaller than 1000x1000. Divide by Plancks constant to get energy in Hz
    [d,sortIdx] = sort(real(D));
    zeemanMap(k,:)=d;        
end
% H = Ham.hyperfine.total;
% [V,D] = eig(full(H/C.h),'vector'); %no need to use sparse matrices for matrices smaller than 1000x1000. Divide by Plancks constant to get energy in Hz
% [d,sortIdx] = sort(real(D));
% v = V(:,sortIdx);
%%
figure(23)
plot(Bs,zeemanMap(:,1:32),color=[1,1,1]*0.5)
%%
zeemanMap(1,1:32)