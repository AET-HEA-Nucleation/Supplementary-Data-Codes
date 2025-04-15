function [indGrain] = meas02grains_ForPtMD_YY_v2_PBC(AtomPos,OP,abcLocalRadius,OPth,CellPara, PBCEnforceYN)


Np = size(AtomPos,2);

% For all sites, local neighbor parent index
indsParent = zeros(1,size(AtomPos,2));
for a1 = 1:Np
    currPos = AtomPos(:,a1);
    
    if OP(a1) < OPth
        indsParent(a1) = -1;
    else
        SubDist = AtomPos - repmat(currPos,[1 size(AtomPos,2)]);
        if PBCEnforceYN
            SubDist(1,:) = mod(SubDist(1,:) + CellPara(1)/2, CellPara(1)) - CellPara(1)/2;
            SubDist(2,:) = mod(SubDist(2,:) + CellPara(2)/2, CellPara(2)) - CellPara(2)/2;
            SubDist(3,:) = mod(SubDist(3,:) + CellPara(3)/2, CellPara(3)) - CellPara(3)/2;
        end
    
        Dist = sqrt( sum( (SubDist).^2, 1));
        subNN = Dist < abcLocalRadius;
        indsNN = find(subNN);
        
        [~,indLocal] = max(OP(subNN));
        indsParent(a1) = indsNN(indLocal);
    end

end

indGrain = zeros(Np,5);
indGrain(:,1:3) = AtomPos';
indGrain(:,5) = OP';

% Traverse grain tree for all sites
for a1 = 1:Np
    if OP(a1) < OPth
        indGrain(a1,4) = -1;
    else
        ind = indsParent(a1);
        search = true;
        while search == true
            if ind == indsParent(ind)
                search = false;
            else
                ind = indsParent(ind);
            end

            if ind==-1
                keyboard
            end
        end
        indGrain(a1,4) = ind;
    end

end

    
end



%end