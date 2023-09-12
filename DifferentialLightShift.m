clear all
addpath("Molecules\");addpath("Atoms\")
C = Constants;
%% Define Fields
maxN = 2;
B.value = 181.5*1e-4; %549.%
B.dir = [0,0,1];
E.value = 0;
E.dir = [0,0,1];
I.value = [0:0.1:6]*1e7;
I.dir = [0,0,1];
I.pol = [1,0,0];
Fields.B = B;
Fields.E = E;
Fields.I = I;
%% Define Hamiltonian
Mol = RbCs;
Ham = Hamiltonian(Molecule = Mol, Fields=Fields, maxN=maxN);
H0 = Ham.hyperfine.total + B.value*Ham.zeeman + E.value*Ham.dc_stark;
H_ac = Ham.ac_stark;
Base = Ham.Basis;
nStates = Base.NStates;
%% Diagonalise for varying intensity
energyMap = nan(length(Fields.I.value),nStates);
statesMap = nan(length(Fields.I.value),nStates,nStates);
for k = 1:length(Fields.I.value)
    H = H0 + Fields.I.value(k)*H_ac;
    [V,D] = eig(full(H/C.h),'vector'); %no need to use sparse matrices for matrices smaller than 1000x1000. Divide by Plancks constant to get energy in Hz
    [d,sortIdx] = sort(real(D));
    energyMap(k,:)=d;        
    statesMap(k,:,:) = V(:,sortIdx);
end
%% Plotting
figure(99)
x = Fields.I.value;
N = Base.getStates("N");
selectStates = N==1;
GSenergy = energyMap(:,1);
offset = 983;
diffEnergy = (energyMap(:,selectStates) - GSenergy - offset*1e6)/1e6;
plot(x,diffEnergy, color=[1,1,1]*0.5)
xlabel("I (W/m^2)")
ylabel(sprintf("E - %dMHz (MHz)",offset))

%% Diagonalise for varying angle
mN = Base.getStates("mN");
Fields.I.value = 2.35e6;
Fields.I.dir = [0,1,0];
angles =  linspace(0,pi/2,50);
Base = Ham.Basis;
nStates = Base.NStates;
energyMapAngle = nan(length(angles),nStates);
statesMapAngle = nan(length(angles),nStates,nStates);
for k = 1:length(angles)
    theta = angles(k);
    Fields.I.pol = [sin(theta),0,cos(theta)];
    H_ac = Ham.makeACStark(Mol.a0,Mol.a2,N,mN,Field=Fields.I);
    H = H0 + Fields.I.value*H_ac;
    [V,D] = eig(full(H/C.h),'vector'); %no need to use sparse matrices for matrices smaller than 1000x1000. Divide by Plancks constant to get energy in Hz
    [d,sortIdx] = sort(real(D));
    energyMapAngle(k,:)=d;        
    statesMapAngle(k,:,:) = V(:,sortIdx);
end
%% Plotting
figure(98)
x = angles*180/pi;
N = Base.getStates("N");
selectStates = N==1;
GSenergy = energyMapAngle(:,1);
offset = 983;
diffEnergyAngles = (energyMapAngle(:,selectStates) - GSenergy - offset*1e6)/1e6;
plot(x,diffEnergyAngles, color=[1,1,1]*0.5)
xlabel("Angle (°)")
ylabel(sprintf("E - %dMHz (MHz)",offset))

%% Reproduce fig 2: robust
B.value = 154.5*1e-4; %549.%
E.value = 0;
Ham = Hamiltonian(Molecule = RbCs, Fields=Fields, maxN=maxN);
h_tmp = Ham.hyperfine;
hyperfine = h_tmp.electricQuadrupole+h_tmp.rigidRotor + h_tmp.spinSpinScalar+ h_tmp.spinRotation;
H0 = hyperfine + B.value*Ham.zeeman + E.value*Ham.dc_stark;
mN = Base.getStates("mN");
I.value = [0:40]*1e7;
I.dir = [0,1,0];
angles =  [0,55,90]*pi/180;
Base = Ham.Basis;
nStates = Base.NStates;
energyMap = nan(length(angles),length(I.value),nStates);
statesMap = nan(length(angles),length(I.value),nStates,nStates);
a0 = 0; %1550nm
a2 = [545]*4*pi*C.e0*C.a0^3; %at 1550nm
for k = 1:length(angles)
    for l = 1:length(I.value)
    theta = angles(k);
    I.pol = [sin(theta),0,cos(theta)];
    H_ac = Ham.makeACStark(a0,a2,N,mN,Field=I);
    H = H0 + I.value(l)*H_ac;
    [V,D] = eig(full(H/C.h),'vector'); %no need to use sparse matrices for matrices smaller than 1000x1000. Divide by Plancks constant to get energy in Hz
    [d,sortIdx] = sort(real(D));
    energyMap(k,l,:)=d;        
    statesMap(k,l,:,:) = V(:,sortIdx);
    end
end
%% Find right states and plot
% selectStates = N==1;
psi0 = Base.prepStates(N=0,mN=0,i1=1.5, mi1 =0.5, i2=3.5,mi2=3.5); %(0,4)_1
psi1 = Base.prepStates(N=0,mN=0,i1=1.5, mi1 =1.5, i2=3.5,mi2=1.5); %(0,3)_0
states = squeeze(statesMap(1,1,:,:));
overlap0 = nan(nStates,1);
overlap1 = nan(nStates,1);
for n = 1:nStates
     stateNow = abs(states(:,n));
     overlap0(n) = stateNow'*psi0;
     overlap1(n) = stateNow'*psi1;
end
[~,idx0] = max(overlap0);
[~,idx1] = max(overlap1);
%%
figure(97);clf;
hold on;
selectStates = N==1;
GSenergy = energyMap(:,:,idx0);
offset = -76e3;

diffEnergy2 = abs(squeeze(energyMap(:,:,idx1) - GSenergy - offset));
plot(I.value,(diffEnergy2(1,:)))
plot(I.value,(diffEnergy2(2,:)))
plot(I.value,(diffEnergy2(3,:)))
