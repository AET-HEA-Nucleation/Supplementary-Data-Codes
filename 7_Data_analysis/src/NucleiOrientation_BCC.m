function [uCorner]=NucleiOrientation(model,blen,lb,ub,dang)
%%
% model=importdata('BOP2models.mat');
% model=model.modelarmps;
% model=model(ico_tid,:)';

%for BCC. lattice: 3*2/sqrt(3)
% lb=3.2; % lower bond length
% ub=3.8; % upper bond length
% blen=3.5; % bond length

%for FCC. lattice: 3*sqrt(2)
% lb=3.4; % lower bond length; FCC.mat
% ub=4.2; % upper bond length
% blen=3.8; % bond length

% dang=5; % angle tolerance from perfect lattice (degree)

del=cosd(90-dang);
Corner=[];
Ncorner=0;
for n=1:size(model,2)
    dif=model-repmat(model(:,n),[1 size(model,2)]);
    dis=sqrt(sum(dif.^2,1));
    dis(n)=1e4;
    ind=find(dis>=lb & dis<=ub);
    if isempty(ind)
        continue;
    end
    for i=1:length(ind)-1
        vec1=model(:,ind(i))-model(:,n);
        vec1=vec1/norm(vec1,2);
        for j=i+1:length(ind)
            vec2=model(:,ind(j))-model(:,n);
            vec2=vec2/norm(vec2,2);
            dotproduct=sum(vec1.*vec2);
            if abs(dotproduct)<=del
                Ncorner=Ncorner+1;
                Corner(Ncorner,:)=[n,ind(i),ind(j)];
            end
        end
    end
end

for n=1:size(Corner,1)
    vec1=model(:,Corner(n,1))-model(:,Corner(n,2));
    vec1=vec1/norm(vec1,2);
    vec2=model(:,Corner(n,1))-model(:,Corner(n,3));
    vec2=vec2/norm(vec2,2);
    vec3=cross(vec1,vec2);
    vec3=vec3/norm(vec3,2);
%     dref=[vec1,-vec1,vec2,-vec2,vec3,-vec3]*blen; %for SC

    dref=[vec1,-vec1,vec2,-vec2,vec3,-vec3,(vec1+vec2+vec3)/2,(vec1+vec2-vec3)/2,(vec1-vec2+vec3)/2,...
        (-vec1+vec2+vec3)/2,(-vec1-vec2+vec3)/2,(-vec1+vec2-vec3)/2,(vec1-vec2-vec3)/2,-(vec1+vec2+vec3)/2]*blen; %for BCC

%     dref=[vec1,-vec1,vec2,-vec2,vec3,-vec3,(vec1+vec2)/2,(vec1-vec2)/2,(-vec1+vec2)/2,(-vec1-vec2)/2,...
%         (-vec2+vec3)/2,(vec2-vec3)/2,(-vec2-vec3)/2,(vec2+vec3)/2,...
%         (vec1+vec3)/2,(vec1-vec3)/2,(-vec1+vec3)/2,(-vec1-vec3)/2]*blen; %for FCC
 
    count=3;
    i=1;
    while i<=count
        ref=dref+repmat(model(:,Corner(n,i)),[1 size(dref,2)]);
        for j=1:size(dref,2)
            dif=model-repmat(ref(:,j),[1 size(model,2)]);
            dis=sqrt(sum(dif.^2,1));
            ind=find(dis<=blen*tand(dang));
            if length(ind)>1
                [~,ind]=min(dis);
            end
            if ~ismember(ind,Corner(n,:))
                count=count+1;
                Corner(n,count)=ind;
            end
        end
        i=i+1;
    end
end
%% identify the locally largest MRO
if isempty(Corner)
    uCorner=[];
    return
end
uind=unique(Corner(:,1));
uCorner=[];
for i=uind'
    ind=find(Corner(:,1)==i);
    nel=sum(Corner(ind,:)~=0,2);
    temp=max(nel);
    temp=nel==temp;
    uCorner=[uCorner;Corner(ind(temp),:)];
end
temp=sum(uCorner~=0,2);
ind=temp>=3; % must have more than 3 atoms
uCorner=uCorner(ind,:);

for i=1:size(uCorner,1)-1
    if uCorner(i,1)==0
        continue
    end
    for j=i+1:size(uCorner,1)
        if uCorner(j,1)==0
            continue
        end
        temp=intersect(uCorner(i,:),uCorner(j,:));
        if length(temp)>=1 %MRO with >=1 common atoms are considered same MRO, and only keep the larger one
            if sum(uCorner(i,:)~=0,2)>sum(uCorner(j,:)~=0,2)
                uCorner(j,:)=0;
            else
                uCorner(i,:)=0;
                break;
            end
            
        end
        
    end
end
uCorner=uCorner(uCorner(:,1)~=0,:);
