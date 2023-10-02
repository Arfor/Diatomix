%% Example - AC Stark Map
clear all
addpath("Molecules\");addpath("Atoms\")
C = Constants;
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot,'defaultAxesBox','on')
set(groot,'defaultAxesLineWidth',0.5)
%% Define Fields
maxN = 2;
B.value = 1;
B.dir = [0,0,1];
B.scaling = 1e-4; %put B.value in Gauss, not Tesla
E.value = 0; %V/m
E.dir = [0,0,1];
E.scaling = 1e2;%put E.value in V/cm, not V/m

Fields.B = B;
Fields.E = E;
%% Define Hamiltonian
Mol = KRb;
I.value = linspace(0,10,51);
I.dir = [0,1,0];
theta=0;
I.pol = [sin(theta),0,cos(theta)];
I.scaling = 1e7;%put I.value in kW/cm^2, not W/m^2
Fields.I = I;
Ham = Hamiltonian(Molecule = Mol, Fields=Fields, maxN=maxN);
H0 = Ham.hyperfine.total + E.value*E.scaling*Ham.dc_stark + B.value*B.scaling*Ham.zeeman;
Base = Ham.Basis;
nStates = Base.NStates;
st = Base.getStates("all");
%% Check matrix elements
% Ham.ac_stark = Ham.makeACStark(0,Mol.a2,st.N,st.mN,"Field",Fields.I);
{min((Ham.ac_stark),[],'all')
max((Ham.ac_stark),[],'all')}
nnz(round(Ham.ac_stark,50)) %get rid of some neglible numbers, left from rotation
% map = [1 0,0;
%     1,1,1;
%     0,0,1];
% figure(3123);
% imagesc(round(Ham.ac_stark,50)); colormap(map)
% axis square
% set(gca,"YDir","normal")
%% Diagonalise for varying intensity
F = I;
x = F.value;
energyMap = nan(length(x),nStates);
statesMap = nan(length(x),nStates,nStates);
disp("Diagonalising - ")
for k = 1:length(x)
    if mod(k,10)==1; fprintf("%d - ",k); end
    H = H0 + x(k)*F.scaling*Ham.ac_stark; %searching for the factor (sqrt(6)/2), either sphericalTensorDot or coupleCartesianSpherically
    [V,D] = eig(full(H/C.h),'vector'); %no need to use sparse matrices for matrices smaller than 1000x1000. Divide by Plancks constant to get energy in Hz
    [d,sortIdx] = sort(real(D));
    energyMap(k,:)=d;        
    statesMap(k,:,:) = V(:,sortIdx);
end
disp("Finished")
%% Plotting
GSEnergy = min(energyMap,[],'all')*1e-6;
GSEnergy = 0;

N = Base.getStates("N");
figure(12); clf;
t = tiledlayout("flow",TileSpacing="tight");

nexttile(t); %N=2 manifold
plot(x,energyMap(:,N==2)*1e-6-GSEnergy, color=[1,1,1]*0.3);
ylabel(sprintf("E (MHz)"))
% xticklabels([])
% ylim([1,4]+5880)
xlim([min(x),max(x)])
axis square
title("N=2")

nexttile(t); %N=1 manifold
plot(x,energyMap(:,N==1)*1e-6-GSEnergy, color=[1,1,1]*0.3);
ylabel(sprintf("E (MHz)"))
% xticklabels([])
% ylim([980,981])
% ylim([3478.8,3480])
xlim([min(x),max(x)])
% xlim([0,3])
axis square
title("N=1")
% 
% nexttile(t); %N=0 manifold
% plot(x,energyMap(:,N==0)*1e-6-GSEnergy, color=[1,1,1]*0.3); hold on;
% ylabel(sprintf("E (MHz)"))
% % ylim([-3,1])
% xlim([min(x),max(x)])
% axis square
% title("N=0")

xlabel("I (kW/cm^2)")
title(t,sprintf("AC Stark Map %s",Mol.name))
subtitle(t, sprintf("E=%.1f V/cm, B=%.1f G", Fields.E.value,Fields.B.value))


