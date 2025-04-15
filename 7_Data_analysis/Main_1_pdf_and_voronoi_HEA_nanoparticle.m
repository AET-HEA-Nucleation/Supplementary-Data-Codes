%% Main_1_pdf_voronoi_BOO_HEA_nanoparticle
% calculate the pair distribution functions and Voronoi analysis
% for the HEA nanoparticles
% This script also provides the option to plot pdf and voronoi using 
% pre-calculated results.
clear;clc;close all
addpath('./src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - calculate pdf and voronoi using 3D atomic model
% 2 - plot pdf and voronoi using pre-calculated results
actionFlag = 2; % Choose 1 or 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if actionFlag == 1
for sampleID = 1:25
%% calculate PDF of all HEA nanoparticles
path='../6_Final_coordinates/';
% read in files: finalized atomic coordinates in Angstrom
data = importdata([path, 'Final_atomic_model_HEA_' num2str(sampleID) '_nanoparticle.mat']);
model = data.model;
% set the parameters: pixel size, step size and the range for pdf
if sampleID==11 || sampleID==15
    res = 0.469;
else
    res = 0.347;
end

step=0.1/res;
cutoff=20/res;

% number of cores to use for parallel computing
parpool_size=12;

model=model/res; % convert to unit of pixels

% calculate the alpha shape of the nanoparticle
submodel=model-repmat(min(model,[],2),[1,size(model,2)])+ones(size(model));
submodel=double(submodel);
shp = alphaShape(submodel(1,:)',submodel(2,:)',submodel(3,:)',ceil(4/res));

% initialize a spherical index for later intersection calculation
[vMat,~]=spheretri(2000);
nsample=size(vMat,1);

% perform the main pdf calculation
xx=0:step:cutoff;
gg=xx*0;
nn=gg;
if parpool_size~=0
pjob = gcp('nocreate');
if isempty(pjob)
    parpool(parpool_size)
elseif pjob.NumWorkers ~= parpool_size
    delete(pjob)
    parpool(parpool_size)
end
end

parfor i=1:size(submodel,2)
    i
    g=xx*0;
    n=g;
    for j=1:size(submodel,2)
        dis=submodel(:,i)-submodel(:,j);
        dis=norm(dis,2);
        if dis<cutoff
            ind=ceil(dis/step+0.01);
            g(ind)=g(ind)+1;
        end
    end
    for j=xx
        spoints=vMat*(j+step/2)+repmat(submodel(:,i)',[size(vMat,1),1]);
        in = inShape(shp,spoints(:,1),spoints(:,2),spoints(:,3));
        ind=round(j/step)+1;
        n(ind)=sum(in)/nsample*(j+step/2)^2;
    end
    gg=gg+g;
    nn=nn+n;
end
rdfnorm=gg./nn;
rdfnorm=rdfnorm/rdfnorm(end-2);
g_r=imgaussfilt(rdfnorm(2:end-1),1.5);
r=res*xx(2:end-1)+0.05;
figure(1);hold on;
plot(r,g_r,'-','linewidth',2)
save(['./Output/PDF_HEA_' num2str(sampleID) '_nanoparticle.mat'],'r','g_r')
figure(1);hold off
%% perform Voronoi analysis
model = data.model';
amorInd=1:size(model,1);    % analyze all atoms

% Voronoi regulations
areacutoff=0.01; % facets area should be larger than 1% of polygon total surface area
valleyposition=importdata('./Input/valleyposition.mat'); % load 1st valley positions of pdf
bondlengthcutoff=valleyposition(sampleID); % Voronoi neighbor should be within the 1st nerest shell distance (1st valley in pdf)

% perform the voronoi calcualtion
defaultFaceColor  = [0.6875 0.8750 0.8984];
[V,R] = voronoin(model);

DT = delaunayTriangulation(model(:,1),model(:,2),model(:,3));
ConnectivityList=DT.ConnectivityList;

vor=[]; vorid=[];   vornk=[];   vorarea=[]; neigh=[];
for tid=1:size(model,1)
tid
XR10 = V(R{tid},:);
if sum(isinf(XR10),[1 2 3])
    ind=sum(ConnectivityList==tid,2)~=0;
    ind=ConnectivityList(ind,:);
    ind=unique(ind(:));
    dis=vecnorm(model(ind,:)-model(tid,:),2,2);
    ind=ind(dis<bondlengthcutoff);
    neigh{tid}=model(ind,:);
    continue
end
K = convhull(XR10);

vec1=XR10(K(:,1),:)-XR10(K(:,2),:);
vec2=XR10(K(:,1),:)-XR10(K(:,3),:);
nk=cross(vec1,vec2);
nk=nk./repmat(sqrt(sum(nk.^2,2)),[1 3]);

ind=1:size(nk,1);
num=0;
faceid=[];
facevor=[];
while ~isempty(ind)
    flag=sum(nk(ind,:).*nk(ind(1),:),2);
    faceidtemp = abs(flag)>1-1e-5 & abs(sum((XR10(K(ind,1),:)-XR10(K(ind(1),1),:)).*nk(ind(1),:),2))<1e-5;
    tempid=K(ind(faceidtemp),:);
    num=num+1;
    faceid{num}=unique(tempid(:));
    
    % sort vertices of each face to follow clockwise or counterclockwise order
    % for plotting
    center=mean(XR10(faceid{num},:),1);
    pol=XR10(faceid{num},:)-center;
    pol=pol./repmat(sqrt(sum(pol.^2,2)),[1 3]);
    npol=size(pol,1);
    Y=dot(cross(repmat(pol(1,:),[npol,1]),pol)',repmat(nk(ind(1),:),[npol,1])');
    
    D = atan2d(Y,dot(repmat(pol(1,:),[npol,1])',pol'));
    D(D<0)=360+D(D<0);
    [~,sid]=sort(D);
    faceid{num}=faceid{num}(sid);
    ind(faceidtemp)=[];
end
facenk=[];
facearea=[];
for i=1:size(faceid,2)
    % calculate surface normal
    vec1=XR10(faceid{i}(1),:)-XR10(faceid{i}(2),:);
    vec2=XR10(faceid{i}(1),:)-XR10(faceid{i}(3),:);
    nk=cross(vec1,vec2);
    facenk(i,:)=nk/sqrt(sum(nk.^2));
    
    % calculate face area
    vec=XR10(faceid{i}(2:end),:)-XR10(faceid{i}(1),:);
    facearea(i)=0.5*sum(vecnorm(cross(vec(1:end-1,:),vec(2:end,:)),2,2));
    facevor{i}=XR10(faceid{i},:);
end
vorid{tid}=faceid;
vor{tid}=facevor;
vornk{tid}=facenk;
vorarea{tid}=facearea;

% remove face with small area
faceremove = find(facearea < areacutoff*sum(facearea));
ind=sum(ConnectivityList==tid,2)~=0;
ind=ConnectivityList(ind,:);
ind=unique(ind(:));
for i=faceremove
    vec=XR10(faceid{i}(1),:)-model(tid,:);
    vec=sum(vec.*facenk(i,:))*facenk(i,:)*2;
    pos=[model(tid,:)+vec;model(tid,:)-vec];
    dis1=vecnorm(model(ind,:)-pos(1,:),2,2);
    dis2=vecnorm(model(ind,:)-pos(2,:),2,2);
    atomremove=find(dis1<1e-5 | dis2<1e-5);
    ind(atomremove)=[];
end
dis=vecnorm(model(ind,:)-model(tid,:),2,2);
% ind=ind(dis<bondlengthcutoff); % comment this line if don't want bond length regulation to the voronoi
neigh{tid}=model(ind,:);
%%%%%%%%%%%%%%% re-calculate voronoi after regulation %%%%%%%%%%%%%%
if length(ind)<5
    vor{tid}=[];
    continue
end
tidsub=find(ind==tid);
[Vsub,Rsub] = voronoin(model(ind,:));
XR10sub = Vsub(Rsub{tidsub},:);
if sum(isinf(XR10sub),[1 2 3])
    vor{tid}=[];
    continue
end
Ksub = convhull(XR10sub);

vec1=XR10sub(Ksub(:,1),:)-XR10sub(Ksub(:,2),:);
vec2=XR10sub(Ksub(:,1),:)-XR10sub(Ksub(:,3),:);
nk=cross(vec1,vec2);
nk=nk./repmat(sqrt(sum(nk.^2,2)),[1 3]);

indsub=1:size(nk,1);
num=0;
faceid=[];
facevor=[];
while ~isempty(indsub)
    flag=sum(nk(indsub,:).*nk(indsub(1),:),2);
    faceidtemp = abs(flag)>1-1e-5 & abs(sum((XR10sub(Ksub(indsub,1),:)-XR10sub(Ksub(indsub(1),1),:)).*nk(indsub(1),:),2))<1e-5;
    tempid=Ksub(indsub(faceidtemp),:);
    num=num+1;
    faceid{num}=unique(tempid(:));
    
    % sort vertices of each face to follow clockwise or counterclockwise order
    % for plotting
    center=mean(XR10sub(faceid{num},:),1);
    pol=XR10sub(faceid{num},:)-center;
    pol=pol./repmat(sqrt(sum(pol.^2,2)),[1 3]);
    npol=size(pol,1);
    Y=dot(cross(repmat(pol(1,:),[npol,1]),pol)',repmat(nk(indsub(1),:),[npol,1])');
    
    D = atan2d(Y,dot(repmat(pol(1,:),[npol,1])',pol'));
    D(D<0)=360+D(D<0);
    [~,sid]=sort(D);
    faceid{num}=faceid{num}(sid);
    indsub(faceidtemp)=[];
    facevor{num}=XR10sub(faceid{num},:);
end
vor{tid}=facevor;
neigh{tid}=model(ind,:);
end
% calculate indices
edgelengthcutoff=0.00; % edge length must be larger than this value

indlist=zeros(1,6); % voronoi index list
Nindlist=0; % counts of voronoi index
VoronoiID=[];
VoronoiID{1}=0; % voronoi cell IDs in each index
badID=[];
facecounts=0; % every voronoi cell has the same weight, not every face
for tid=amorInd % only count amorphous region
    facevor=vor{tid};
    if isempty(facevor)
        continue
    end
    edgelength=[];
    nedge=0;
    totedgelen=0;
    for i=1:size(facevor,2)
        ntemp=size(facevor{i},1);
        nedge=nedge+ntemp;
        edgelength{i}=vecnorm(facevor{i}-facevor{i}([ntemp 1:ntemp-1],:),2,2);
        if sum(edgelength{i}>6)
            badID=[badID tid];
            continue
        end
        totedgelen=totedgelen+sum(edgelength{i});
    end
    ind=zeros(1,6);
    for i=1:size(facevor,2)
        n=sum(edgelength{i}>=totedgelen/nedge*edgelengthcutoff);
        if n<=6
            ind(n)=ind(n)+1;
        end
    end
    facecounts=facecounts+ind/size(neigh{tid},1);
    
    temp=indlist==ind;
    id=find(sum(temp,2)==6);
    if ~isempty(id)
        Nindlist(id)=Nindlist(id)+1;
        VoronoiID{id}=[VoronoiID{id} tid];
    else
        indlist=[indlist;ind];
        Nindlist=[Nindlist,1];
        VoronoiID{end+1}=[tid];
    end
end
save(['./Output/Voronoi_HEA_' num2str(sampleID) '_nanoparticle.mat'],'Nindlist','indlist','neigh')
end
else % plot pdf and voronoi using pre-calculated results
%% plot PDF of HEA nanoparticles
peakposition = importdata(['./Input/peakposition.mat']);
nind=0;
figure(1)
clf
for sampleID=1:25
nind=nind+1;
figure(1)
hold on
data=importdata(['./Input/PDF_HEA_' num2str(sampleID) '_nanoparticle.mat']);

xx=data.r/peakposition(sampleID);
g_r=data.g_r;

c1 = [0,0.447,0.741];
c2 = [0.85,0.325,0.098];

plot(xx,g_r+nind*0.5,'linewidth',1,'color',c1+(c2-c1)*nind/26);
end

ylim([0,18]);xlim([0 6.5])
hold on;ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1)
hold off
xlabel(['r/R_1'],'fontsize',14)
ylabel('g(r)','fontsize',14)
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold')
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;

%% plot Voronoi of HEA2，6，10 nanoparticles
pltlt = [2 6 10]; % select datasets to plot

data=[];ind=[];NNico=[];N=[];nfold = [];nfoldplt = [];
for i=1:length(pltlt)
    data{i}=importdata(['./Input/Voronoi_HEA_' num2str(pltlt(i)) '_nanoparticle.mat']);
    ind{i}=data{i}.indlist;
    N{i}=data{i}.Nindlist;

    nfold{i}=ind{i}.*N{i}';
    nfoldplt(:,i)=sum(nfold{i}(:,3:6),1)/sum(nfold{i}(:))';

    idico=ind{i}(:,5)>=8 & abs(sum(ind{i},2)-12)<=1;
    NNico(1,i)=sum(N{i}(idico))/sum(N{i});
end

NNN=12; % number of voronoi indices to plot

pltN = 1;% order index based on the first particle
[~,indt]=sort(N{pltN},'descend');
reftop10=ind{pltN}(indt(1:NNN),:);

% find counts in particles
Ncount=[];NN = [];
for i=1:NNN
    for j=1:length(pltlt)
        id=find(sum(ind{j}==reftop10(i,:),2)==6);
        if isempty(id)
            Ncount{j}(i)=0;
        else
            Ncount{j}(i)=N{j}(id);
        end
    end
end
for j=1:length(pltlt)
    Ncount{j}=Ncount{j}/sum(N{j});
    NN=[NN Ncount{j}'];
end

% order based on coordination number
[~,indt]=sort(N{pltN},'descend');
top10=ind{pltN}(indt(1:NNN),3:end);
topOrder=sum(top10,2)*1e5+top10*([1e4 1e3 1e2 1e1]');
[~,indt]=sort(topOrder,'ascend');

figure(2) % plot voronoi diagram for selected particles
b=bar(NN(indt(1:end),:)*100,1);
b(1).FaceColor=[0 0.5 0];
b(2).FaceColor=[0, 0.4470, 0.7410];
b(3).FaceColor=[0.6350, 0.0780, 0.1840];
box off
ylabel('Fraction (%)')
temp=['<'*ones(NNN,1) num2str(top10(indt,:)) '>'*ones(NNN,1)];
label=cell(NNN,1);
for i=1:NNN
    label{i}=temp(i,:);
end
label = cellfun(@(x) strrep(x,'  ',','), label,'UniformOutput',false);
set(gca,'xticklabel',label);xtickangle(90);ax = axis;hold on
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1);hold off
legend(['HEA' num2str(pltlt(1))], ['HEA' num2str(pltlt(2))],['HEA' num2str(pltlt(3))]) % Add a legend
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold')
legend boxoff;ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;
end
