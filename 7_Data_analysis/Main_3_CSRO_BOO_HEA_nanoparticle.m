%% Main_3_CSRO_BOO_HEA_nanoparticle
% Calculate the chemical short range order parameters for the HEA nanoparticles
% and their correlation with crystallinity.
% This script also provides the option to plot CSRO and their BOO using 
% pre-calculated results.
clear;clc;close all
addpath('./src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - calculate chemical order of HEA nanoparticles
% 2 - plot correlation between chemical order and crystallinity using pre-calculated results
actionFlag = 2; % Choose 1 or 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if actionFlag == 1
for sampleID = 1:25
%% calculate chemical short range order parameters of all HEA nanoparticles

% read in files: finalized atomic coordinates in Angstrom and BOO 
data = importdata(['./Input/BOO_HEA_' num2str(sampleID) '_nanoparticle.mat']);
order=data.order;
model=data.model;
atom=data.atom;
% set the parameters: pixel size
if sampleID==11 || sampleID==15
    res = 0.469;
else
    res = 0.347;
end
% define neighbor using voronoi 
dataVor=importdata(['./Input/Voronoi_HEA_' num2str(sampleID) '_nanoparticle.mat']);
neigh=dataVor.neigh;
neighind=cell(size(model,2),1);
for i=1:size(neigh,2)
    for j=1:size(neigh{i},1)
        dif=sum(abs(model-repmat(neigh{i}(j,:)',[1 size(model,2)])),1);
        neighind{i}=[neighind{i}; find(dif<1e-2)];
    end
    ind=neighind{i};
    ind(ind==i)=[];
    neighind{i}=ind;
end

% calculate global chemical order
C=zeros(3,1);
for i=1:3
    C(i)=sum(atom==i)/length(atom);
end
alpha1=zeros(3,3);
for i=1:3
    indi=find(atom==i);
    ind=[];
    for k=indi
        ind=[ind;neighind{k}];
    end
    atomt=atom(ind);
    for j=1:3
        alpha1(i,j)=(sum(atomt==j)/length(ind)-C(j))/((i==j)-C(j));
%         alpha1(i,j)=(sum(atomt==j)/sum(atom==i)/size(nei_vec,2)-C(j))/((i==j)-C(j));
%         alpha1(i,j)=(sum(atomt==j)/sum(atom==i)/12-C(j))/((i==j)-C(j));
    end
end
alpha=(alpha1+permute(alpha1,[2 1 3]))/2;
% calculate CSRO for each atom
radius=5;
alpha_arr=zeros(3,3,length(atom));
C=zeros(3,1);
for i=1:3
    C(i)=sum(atom==i)/length(atom);
end

for n=1:length(atom)
    dis=model-repmat(model(:,n),[1 length(atom)]);
    dis=sqrt(sum(dis.^2,1));
    atom_ind_radius=dis<radius;
for i=1:3
    indi=find(atom==i & atom_ind_radius);
    ind=[];
    for k=indi
        ind=[ind;neighind{k}];
    end
    atomt=atom(ind);
    for j=1:3
%         alpha_arr(i,j,n)=(sum(atomt==j)/sum(atom==i & atom_ind_radius)/12-C(j))/((i==j)-C(j));
        alpha_arr(i,j,n)=(sum(atomt==j)/length(ind)-C(j))/((i==j)-C(j));
    end
end
end
alpha_arr(isnan(alpha_arr))=0;
alpha_arr=(alpha_arr+permute(alpha_arr,[2,1,3]))/2;
save(['./Output/CSRO_HEA_' num2str(sampleID) '_nanoparticle.mat'], 'alpha_arr','alpha1','alpha')
end
else % plot correlation between chemical order and crystallinity using pre-calculated results
%% plot correlation between chemical order and crystallinity
% load CSRO and BOO results
order_tot=[];alpha_arr_tot=[];
for i=1:25
    data=importdata(['./Input/CSRO_HEA_' num2str(i) '_nanoparticle.mat']);
    alpha_arr = data.alpha_arr;
    data2=importdata(['./Input/BOO_HEA_' num2str(i) '_nanoparticle.mat']);
    order=data2.order;
    model=data2.model;
    atom=data2.atom;
    order_tot = cat(2,order_tot,order-mean(order));
    alpha_arr_tot = cat(3,alpha_arr_tot,alpha_arr);
end
order = order_tot;
alpha_arr = alpha_arr_tot;
% remove sigularity values
temp=squeeze(alpha_arr(1,1,:))';
ind=temp==0;alpha_arr(:,:,ind)=[];order(ind)=[];
temp=squeeze(alpha_arr(2,2,:))';
ind=temp==0;alpha_arr(:,:,ind)=[];order(ind)=[];
temp=squeeze(alpha_arr(3,3,:))';
ind=temp==0;alpha_arr(:,:,ind)=[];order(ind)=[];
% calculate average BOO values for different ranges of CSRO
del = 0.2;x_range = -1:del:1;
BOO_avg = zeros(6,length(x_range));
i=1;
for x = x_range
    ind = abs(alpha_arr-x)<=del;
    BOO_avg(1,i)=mean(order(ind(1,1,:)));
    BOO_avg(2,i)=mean(order(ind(2,2,:)));
    BOO_avg(3,i)=mean(order(ind(3,3,:)));
    BOO_avg(4,i)=mean(order(ind(2,3,:)));
    BOO_avg(5,i)=mean(order(ind(1,3,:)));
    BOO_avg(6,i)=mean(order(ind(1,2,:)));
    i=i+1;
end
% alpha_33 vs BOO
figure(1)
bar(x_range,BOO_avg(3,:))
ylabel('\DeltaBOO');xlabel('CSRO')
xlim([-1.1 1.1]);ylim([-0.15 0.15]);legend('\alpha_{33}')
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold');box on
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;
% all six alphas vs BOO
figure(2)
subplot(121)
bar(x_range,BOO_avg(1:3,:));
ylabel('\DeltaBOO');xlabel('CSRO');legend('\alpha_{11}','\alpha_{22}','\alpha_{33}')
xlim([-1 1]*1.1);ylim([-0.15 0.150001])
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold');box on
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;
subplot(122)
bar(x_range,BOO_avg([6 5 4],:))
ylabel('\DeltaBOO');xlabel('CSRO');legend('\alpha_{12}','\alpha_{13}','\alpha_{23}')
xlim([-1 1]*1.1);ylim([-0.15 0.150001])
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold');box on
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;
end
