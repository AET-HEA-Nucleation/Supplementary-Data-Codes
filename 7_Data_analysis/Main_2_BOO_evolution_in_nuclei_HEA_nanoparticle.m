%% Main_2_BOO_evolution_in_nuclei_HEA_nanoparticle
% Calculate the bond orientation order parameters for the HEA nanoparticles
% and identify nuclei as well as their BOO evolutions.
% This script also provides the option to plot nuclei and their BOO using 
% pre-calculated results.
clear;clc;close all
addpath('./src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - calculate BOO and identify nuclei using 3D atomic model
% 2 - plot nuclei and their BOO using pre-calculated results
actionFlag = 2; % Choose 1 or 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if actionFlag == 1
for sampleID = 1:25
%% calculate BOO of all HEA nanoparticles
path='../6_Final_coordinates/';
% read in files: finalized atomic coordinates in Angstrom
data = importdata([path, 'Final_atomic_model_HEA_' num2str(sampleID) '_nanoparticle.mat']);
model = data.model;
atoms = data.atoms;
% set the parameters: pixel size, step size and the range for pdf
if sampleID==11 || sampleID==15
    res = 0.469;
else
    res = 0.347;
end
valleyposition=importdata('./Input/valleyposition.mat'); % load 1st valley positions of pdf
cutoff = valleyposition(sampleID); % first valley position in pdf

% calculte averaged local bond orientation order parameter
Q4bar=qnbar(4,model,cutoff);
Q6bar=qnbar(6,model,cutoff);

% reference values of crystalline lattice
Q4barb=0.0363696; Q6barb=0.510688;
Q4barf=0.190941;Q6barf= 0.574524;
Q4barh=0.09722; Q6barh=0.484762;

% normalize BOO by using fcc lattice
order=sqrt(Q4bar.^2+Q6bar.^2)/sqrt(Q4barf.^2+Q6barf.^2);
save(['./Output/BOO_HEA_' num2str(sampleID) '_nanoparticle.mat'],'model','order','atoms')

%% search crystal nuclei
OPth = 0.5; % BOO threshold for nuclei atoms

data=importdata(['./Output/BOO_HEA_' num2str(sampleID) '_nanoparticle.mat']);
scaled_SROP=data.order;
model_OriOri=double(data.model);
atoms=data.atoms;

% maximum separation between a nucleus atom and its belonging nucleus
abcLocalRadius = 4;   PBCEnforceYN = 0;    CellPara = [];   PerfCase = 1;

% regulations on the BOO singlarities values
Q_tot = scaled_SROP;
OP = Q_tot;  minOP = min(OP);
OP(OP > PerfCase) = PerfCase + (PerfCase - OP(OP > PerfCase));
OP(OP < minOP) = minOP;  OP(isnan(OP)) = minOP;

% assign nuclei atoms based on atoms with locally highest BOO
[indGrain] = meas02grains_ForPtMD_YY_v2_PBC(model_OriOri,OP,abcLocalRadius,OPth,CellPara, PBCEnforceYN);
fprintf('crystal region atom num: %i\n',sum(indGrain(:,4)~=-1))

% apply 13 atom constraint
indGrain_ori=indGrain;
indGrain_type = unique(indGrain(:,4));
indGrain_type(1) = [];
for i = 1:numel(indGrain_type)
    if sum(indGrain(:,4)==indGrain_type(i))<13
        indGrain(indGrain(:,4)==indGrain_type(i),4) = -1;
    end
end
fprintf('crystal nuclei atom num: %i\n',sum(indGrain(:,4)~=-1))

model_nuclei = model_OriOri(:,indGrain(:,4)~=-1);
atoms_nuclei = atoms(indGrain(:,4)~=-1);

% identify nuclei orientation based on brute-force seearch and template matching
peakarray=importdata('./Input/peakposition.mat');
valleyposition=importdata('./Input/valleyposition.mat');

typegrain = unique(indGrain_ori(:,4));
typegrain(1)=[];

indGrain_combined=indGrain_ori;
indGrain_combined(:,6)=-1;

blen=peakarray(sampleID)*sqrt(2); % optimum FCC lattice constant by pdf 1st peak
lb=blen*0.9; % lower bound for lattice constant
ub=blen*1.1; % upper boound for lattice constant
dang=5; % maxium allowed angle torlance (deg)
del=0.15; % max allowed position error from a fitted lattice point, in percentage of FCC lattice parameter
dang2=5; % nulcei orientation mismatch tolorance

count=0;vector=[];FCCind=[];vectorfit=[];FCCindfit=[];bondfit=[];orifit=[];
for i=typegrain'
    ind=find(indGrain_ori(:,4)==i);
    indmodel=model(:,ind);
    [uCorner]=NucleiOrientation(indmodel,blen,lb,ub,dang);
    fprintf('crystal region atom num by tree search: %i\n',length(ind))
    fprintf('identified atom num by orientation vector search: %i\n',size(uCorner,2))
    if size(uCorner,2)==0
        indGrain_combined(ind,6)=-1;
    else
        [~,temp]=max(sum(uCorner~=0,2));
        uCorner=uCorner(temp,:);
        count=count+1;
        indGrain_combined(ind,6)=count;
        
        vec1=indmodel(:,uCorner(1))-indmodel(:,uCorner(2));
        vec1=vec1/norm(vec1,2);
        vec2=indmodel(:,uCorner(1))-indmodel(:,uCorner(3));
        vec2=vec2/norm(vec2,2);
        vec3=cross(vec1,vec2);
        vec3=vec3/norm(vec3,2);
        vector(:,:,count)=[vec1 vec2 vec3];
        FCCind{count}=uCorner;
               
        [uFCCind,vec,blenFit,ori]=NucleiOrientationFit(indmodel,uCorner,blen,del);
        
        vectorfit(:,:,count)=vec;
        FCCindfit{count}=uFCCind;
        bondfit(count)=blenFit;
        orifit(:,count)=ori;
    end
end

% combine nuclei with similar orientation
fprintf('Total crystal region atom num by tree search: %i\n',sum(indGrain_combined(:,4)~=-1))
fprintf('Total crystal region atom num by orientation vector search: %i\n',sum(indGrain_combined(:,6)~=-1))

combine=[];
for i=1:count-1
    for j=i+1:count
        ind1=find(indGrain_combined(:,6)==i);
        ind2=find(indGrain_combined(:,6)==j);
        flg1=false;
        model1=model(:,ind1);
        model2=model(:,ind2);
        
        for k=1:length(ind1)
            dif=model2-repmat(model1(:,k),[1 length(ind2)]);
            dis=sqrt(sum(dif.^2,1));
            if sum(dis<=valleyposition(sampleID))
                flg1=true;
                break
            end
        end
        flg2=true;
        mat1=vectorfit(:,:,i);
        mat2=vectorfit(:,:,j);
        for i1=1:3
            for j1=1:3
                if abs(sum(mat1(:,i1).*mat2(:,j1)))<=cosd(90-dang2) || abs(sum(mat1(:,i1).*mat2(:,j1)))>=cosd(dang2)
                else
                    flg2=false;
                    break
                end
            end
        end
        
        if flg1 && flg2
            n=find(sum(combine==i | combine==j,2));
            if isempty(n)
                temp=zeros(1,size(combine,2));
                temp(1:3)=[0 i j];
                combine=[combine;temp];
            elseif length(n)==1
                temp=unique([combine(n,:) i j]);
                combine(n,1:length(temp))=temp;
            else
                temp=combine(n,:);
                temp=unique([temp(:)' i j]);
                combine(n(1),1:length(temp))=temp;
                combine(n(2:end),:)=[];
            end
        end
    end
end
indGrain_combined(:,7)=indGrain_combined(:,6);
for i=1:size(combine,1)
    temp=unique(combine(i,:));
    temp(1)=[];
    for j=temp
        ind=indGrain_combined(:,6)==j;
        indGrain_combined(ind,7)=temp(1);
    end
end
save(['./Output/Nuclei_atoms_HEA_' num2str(sampleID) '_nanoparticle.mat'], ...
    'indGrain','indGrain_ori','indGrain_combined','vector','FCCind','vectorfit','FCCindfit','bondfit','orifit')
end

%% analyze BOO evolution from nuclei cores to surfaces of different sizes
xx=0:0.5:20; % distances from nuclei cores to be analyzed
ngauss=4;[vMat,~]=spheretri(100); % spatial interpolation parameters

Nparticle = 25; % total number of HEA particles
Nsize = cell(1,Nparticle);

sizeMode = 1; % 1 - number of atoms; 2 - total BOO
size_arr = [1 10; 10 100; 100 500; 500 1000; 1000 1e4];

size_arr_dim = size(size_arr,1);
order_dis_arr = cell(Nparticle,size_arr_dim);
parfor Nii=1:Nparticle
Nii
% load nuclei information and 3D atomic models
data=importdata(['./Input/Nuclei_atoms_HEA_' num2str(Nii) '_nanoparticle.mat']);
indGrain=data.indGrain_ori; 

data=importdata(['./Output/BOO_HEA_' num2str(Nii) '_nanoparticle.mat']);
model = data.model;
order = data.order;
order(order>1) = 2-order(order>1);

% calculate size of each nucleus
uniqueID = unique(indGrain(:,4));uniqueID(1)=[];
Nsize{Nii} = zeros(2,length(uniqueID));
j=0;
for i = uniqueID'
    j=j+1;
    Nsize{Nii}(1,j)=sum(indGrain(:,4)==i);
    Nsize{Nii}(2,j)=sum(indGrain(indGrain(:,4)==i,5));
end
% calculate BOO evolution of nucleus by interpolation
for jj = 1:size_arr_dim
tempIDarr = find(Nsize{Nii}(sizeMode,:)>=size_arr(jj,1) & Nsize{Nii}(sizeMode,:)<=size_arr(jj,2));
order_dis_arr{Nii,jj} = zeros(length(tempIDarr),length(xx));
j=0;
for tempID = tempIDarr
j=j+1;
nucleiInd = uniqueID(tempID);
% center=indGrain(nucleiInd,1:3); % choose the definition of center
nuc = indGrain(indGrain(:,4)==nucleiInd,1:3);center = mean(nuc,1);

dis = model-center';
dis = sqrt(sum(dis.^2,1));
ind = dis<= max(xx)+ngauss+1;
interpOBJ=rbfcreate(model(:,ind), order(ind),'RBFFunction', 'linear', 'RBFConstant', ngauss);
shp = alphaShape(model(1,ind)',model(2,ind)',model(3,ind)',ngauss);

order_dis = zeros(1,length(xx));
for i = 1:length(xx)
    pos_arr=center+vMat*xx(i);
    in = inShape(shp,pos_arr(:,1),pos_arr(:,2),pos_arr(:,3));

    order_vec = rbfinterp(pos_arr(in,:)', interpOBJ);
    if i>1
        ind = order_vec<order_dis(i-1);
        order_dis(i) = mean(order_vec(ind));
    else
        order_dis(i) = mean(order_vec);
    end
end

order_dis(isnan(order_dis))=0;
order_dis_arr{Nii,jj}(j,:)=order_dis;
end
end
end
save(['./Output/Nuclei_statistics_BOOvsSize_HEA_nanoparticle.mat'], 'xx','order_dis_arr','size_arr')

else % plot nuclei and their BOO using pre-calculated results
%% plot nuclei size distribution
SizeStat = [];count = 0;
for sampleID=1:25
data=importdata(['./Input/Nuclei_atoms_HEA_' num2str(sampleID) '_nanoparticle.mat']);
indGrain=data.indGrain_combined;
temp=sort(unique(indGrain(:,4)));
temp(1)=[];
for jj=temp'
count=count+1;
SizeStat(count,:) = [sampleID jj sum(indGrain(:,4)==jj)];
end
end

step = 20;
figure(1)
[x,y]=hist(SizeStat(:,3),step/2:step:1300);
bar(y,x);set(gca,'YScale','log')
ylim([8e-1,1e4]);xlim([0 1300]);xlabel('Nuclei size');ylabel('Number of nuclei')
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold')
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;
%% plot nuclei size and core BOO
BOOSize = [];
count = 0;
for sampleID=1:25
data=importdata(['./Input/Nuclei_atoms_HEA_' num2str(sampleID) '_nanoparticle.mat']);
indGrain=data.indGrain_combined;
temp=sort(unique(indGrain(:,4)));
temp(1)=[];
for jj=temp'
count=count+1;
BOOSize(count,:) = [sampleID jj indGrain(jj,5) sum(indGrain(:,4)==jj)];
end
end

xstep=0.02;x = 0.5:xstep:1;y = zeros(size(x));j=0;
for i=x
    j=j+1;
    ind = abs(BOOSize(:,3)-i)<=xstep/2;
    y(j) = mean(BOOSize(ind,4));
end
figure(2)
plot(x,y,'-','linewidth',2);axis tight on
xlabel('BOO of nuclei core','FontSize',14);ylabel('Avg. nuclei size','FontSize',14)
set(gca,'FontSize',14,'FontName', 'Arial');hold on
plot([0.49 1.01], [0.95 0.95],'k-','linewidth',1.5);hold off
xlim([0.49 1.01]);ylim([1 300]);set(gca,'Yscale','log');grid off;box on
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold')
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;

%% plot BOO evolution from core to surface of nuclei with different sizes and fit with an exp funtion
load('./Input/Nuclei_statistics_BOOvsSize_HEA_nanoparticle.mat')
Nparticle = 25;
% average BOO evolution curves of nuclei with similar sizes
order_dis_mean = zeros(size(size_arr,1),length(xx));
for jj = 1:size(size_arr,1)
    temp=[];
for Nii=1:Nparticle
    temp = [temp;order_dis_arr{Nii,jj};];
end
order_dis_mean(jj,:) = sum(temp,1)./sum(temp~=0,1);
end
% plot BOO evolution and fit with an exp function
abc=[];
figure(3);clf
scatter(xx,order_dis_mean,'o','filled');hold on
for n = 1:size(size_arr,1)
ind = order_dis_mean(n,:)>=0.51;
x0 = xx(ind);y0=order_dis_mean(n,ind);
    fo = fitoptions('Method','NonlinearLeastSquares','algorithm','Levenberg-Marquardt','TolFun',10^-10,'Robust','Bisquare','TolX',10^-6,...
               'Lower',[y0(1),8,1],...
               'Upper',[y0(1),18,3.5],...
               'StartPoint',[y0(1) 13 1.4]);
ft = fittype('a.*exp(-(x./b).^c)','options',fo);
[curve2,gof2] = fit(x0',y0',ft);
xplot=0:0.04:12;
yplot=curve2.a.*exp(-(xplot./curve2.b).^curve2.c);
abc(n,:) = [curve2.a curve2.b curve2.c];
plot(xplot,yplot,'-','color',[1 1 1]*0.5);
end
hold off
xlabel('Distance to nuclei core (A)');ylabel('BOO');
xlim([0 10.9]);ylim([0.5 1]);legend('1-10','10-100','100-500','500-1000','>1000')
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold');box on
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;

% calculate energy vs nuclei sizes based on gradient nucleation pathway and
% fitted BOO curves
figure(4);clf;
for n=1:5
dg=1;gamma=1;N=20;
lambda=0.01:0.01:2.5;
DG=zeros(size(lambda));
i=0;
for lam=lambda
    i=i+1;
    r=(0:lam/N:lam*5)+lam/N/2;
    DG(i) = -dg*4*pi*(r(2)-r(1)).*sum(r.^2.*abc(n,1).*exp(-(r./lam).^abc(n,3)))-...
        gamma*4*pi*sum(r.^2.*(abc(n,1).*exp(-((r+lam/N/2)./lam).^abc(n,3))-abc(n,1).*exp(-((r-lam/N/2)./lam).^abc(n,3))));
end
plot(lambda,DG,'-','linewidth',2);hold on
end
xlabel('R (\gamma/\Deltag)');ylabel('Energy (\gamma^3/\Deltag^2)');
xlim([0 max(lambda)]);ylim([-1 1]*15)
legend('1-10','10-100','100-500','500-1000','>1000')
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold')
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;hold off
end
