function [diabat, diabatStates, adiabat, adiabatStates] = findAdiabat(x,y,states, xIdx,yIdx)
%[adiabat, adiabatStates, diabat, diabatStates] = findAdiabat(x,y,states, xIdx,yIdx)
% x = array of x values
% y = spectrum for certain number of states
% states = array containing the state-composition of every state in y
% xIdx,yIdx = starting indices to look for adiabat 
% adiabats are under the constraint that the mF cannot change
s = size(states);
nStates = s(end);
%Checks adjacent states to determine least change in state
adiabat = y(:,yIdx); %follow xth energy state, stores energies
adiabatStates = abs(squeeze(states(:,:,yIdx))).^2;

diabatYidx = nan(1,length(x)); %stores yIdx, not energy
diabat = nan(1,length(x)); %stores yIdx, not energy

diabatYidx(xIdx) = yIdx;
diabat(xIdx) = y(xIdx,yIdx);
diabatStates = nan(size(adiabatStates));
diabatStates(xIdx,:) = (squeeze(states(xIdx,:,yIdx)));
% [~, sIdx] = sortrows(diabatStates(xIdx,:)'.^2,"descend");
adjC = 3; %set maximum number of adjacent states to look at
for k=xIdx:(length(x)-1) %forward loop
    stateNow = abs(squeeze(states(k,:,diabatYidx(k))));
    stateSearch =diabatYidx(k)-adjC:diabatYidx(k)+adjC; stateSearch(stateSearch<=0)=[];
    statesNext = abs(squeeze(states(k+1,:,stateSearch)));
    [~,S]=sort(sum((statesNext.*stateNow'),1)); %computes overlap integral
    idx = stateSearch(S(end)); %restrict search to adjacent states
    if abs(idx-diabatYidx(k))>3
        idx = stateSearch(S(end-1));
    end
    diabatYidx(k+1) = idx;
    diabat(k+1) = y(k+1,idx);
    diabatStates(k+1,:) = (squeeze(states(k+1,:,idx)));
end
for k=(xIdx-1):-1:1 %backward loop
    stateIdx = diabatYidx(k+1);
    stateNow = abs(squeeze(states(k+1,:,stateIdx)));
    stateSearch =stateIdx-adjC:stateIdx+adjC; 
    stateSearch(stateSearch<=0)=[];stateSearch(stateSearch>nStates)=[];
    statesPrevious = abs(squeeze(states(k,:,stateSearch))); 
    [~,S]=sort(sum((statesPrevious.*stateNow'),1)); %maximum overlap
    idx = stateSearch(S(end));
    if abs(idx-stateIdx)>3
        idx = stateSearch(S(end-1));
    end
    diabatYidx(k) = idx;
    diabat(k) = y(k,idx);
    diabatStates(k,:) = (squeeze(states(k,:,idx)));
end    

end

