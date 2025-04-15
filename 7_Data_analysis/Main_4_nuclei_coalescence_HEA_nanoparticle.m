%% Main_4_nuclei_coalescence_HEA_nanoparticle
% Calculate the misorientation and distance between all nuclei pairs in the
% HEA nanoparticles, and identify different types of twins based on
% coincidence site lattice theory.
% This script also provides the option to plot nuclei coalescence statistics using 
% pre-calculated results.
clear;clc;close all
addpath('./src')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - calculate the misorientation and distance between all nuclei pairs in HEAs
% 2 - plot nuclei coalescence statistics using pre-calculated results
actionFlag = 2; % Choose 1 or 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if actionFlag == 1
%% calculate the misorientation and distance between all nuclei pairs in HEAs
Nsize_cutoff=13; % only count template larger than 13 atoms
nucleiang=cell(1,25);
parfor sampleID = 1:25
count=0;
sampleID

% read in files: nuclei orientation information
data = importdata(['./Input/Nuclei_atoms_HEA_' num2str(sampleID) '_nanoparticle.mat']);
indGrain_combined=data.indGrain_combined;
FCCind=data.FCCindfit;
vector=data.vectorfit;

% calculate misorientation and distance between all nuclei pairs
typeind=unique(indGrain_combined(:,6));typeind(1)=[];
for i=1:length(typeind)-1
    if length(FCCind{i})>Nsize_cutoff
    for j=i+1:length(typeind)
        if length(FCCind{j})>Nsize_cutoff
        m1=indGrain_combined(indGrain_combined(:,6)==i,1:3);
        m2=indGrain_combined(indGrain_combined(:,6)==j,1:3);
        c1=indGrain_combined(indGrain_combined(:,6)==i,5);
        c2=indGrain_combined(indGrain_combined(:,6)==j,5);
        c = mean([c1;c2]);
        mdis=1e4;
        for k=1:size(m1,1)
            dis=m2-m1(k,:);
            dis=sqrt(sum(dis.^2,2));
            mdis=min([dis(:); mdis]);
        end
        v1=vector(:,:,i);
        v2=vector(:,:,j);
        ind = perms(1:3);
        dir = [];
        for ii = 1:size(ind,1)
            for i1 = 0:1
                for i2 = 0:1
                    v2t = v2(:,ind(ii,:))*[(-1)^i1 0 0;0 (-1)^i2 0;0 0 1];
                    if det(v2t)<0
                        v2t(:,3) = -v2t(:,3);
                    end
                    temp=v1*v2t';
                    dir(end+1,:)=rotm2axang(temp)';
                end
            end
        end
        [~,ind] = min(dir(:,4));
        count=count+1;
        nucleiang{sampleID}(count,:)=[sampleID i j mdis [0 0 0]*180/pi dir(ind,1:3) dir(ind,4)/pi*180 length(FCCind{i}) length(FCCind{j}) mean(c1) mean(c2) c];
        end
    end
    end
end
end
save(['./Output/Nuclei_coalescence_HEA_nanoparticle.mat'], 'nucleiang')
else % plot nuclei coalescence statistics using pre-calculated results
%% plot misorientation angle vs pair separation
nucleiangCell = importdata('./Input/Nuclei_coalescence_HEA_nanoparticle.mat');
nucleiang = cell2mat(nucleiangCell'); % convert cell to matrix in matlab

figure(1)
hist3(nucleiang(:,[4 11])./[10 1],'ctrs',{0:0.1:18 1:1:75},'CDataMode','auto','FaceColor','interp');
myColorMap = jet(256*2)*1;myColorMap(1,:)=1;colormap(myColorMap);axis tight on
xlabel('Pair separation (nm)','FontSize',14)
ylabel('\theta (^\circ)','FontSize',14)
zlabel('Number of nuclei pairs','FontSize',14)
set(gca,'FontSize',14,'FontName', 'Arial');hold off
xlim([0 15]);ylim([0 61]);view(-50,30)
grid off;ax = gca;box on
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold')
ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;

%% plot rotation axis stereograph for misorientation angle of 0, 37, and 60 deg
valleyposition=importdata('./Input/valleyposition.mat'); % load 1st valley positions of pdf
dis_cut = mean(valleyposition);
ang_cut_arr = [0 37 60]; % misorientation angle of 0, 37, and 60 deg
ang_rang_cut = 5; % angle tolerance range

dis_arr = [20 0]; %distance array (angstrom)
% quantify misorientation information of selected nuclei pairs
for ijk1= [1 2 3]
    for ijk2 = 1:2 
ang_cut=ang_cut_arr(ijk1);
nu_ang=nucleiang(abs(nucleiang(:,4)-dis_arr(ijk2))<=dis_cut,:);
ind=abs(nu_ang(:,11)-ang_cut)<=ang_rang_cut;
nu_ang=nu_ang(ind,:);
vec_arr = [];
for iii = 1:size(nu_ang,1)
kk=nu_ang(iii,1);i1=nu_ang(iii,2);j1=nu_ang(iii,3);
if iii == 1
data=importdata(['./Input/Nuclei_atoms_HEA_' num2str(kk) '_nanoparticle.mat']);
indGrain_combined=data.indGrain_combined;
FCCind=data.FCCindfit;
vector=data.vector;
elseif kk ~=nu_ang(iii-1,1)
data=importdata(['./Input/Nuclei_atoms_HEA_' num2str(kk) '_nanoparticle.mat']);
indGrain_combined=data.indGrain_combined;
FCCind=data.FCCindfit;
vector=data.vector;
end
vec = nu_ang(iii,8:10)*vector(:,:,i1);
vec1 = sort(abs(vec/norm(vec)),'descend');
vec = nu_ang(iii,8:10)*vector(:,:,j1);
vec2 = sort(abs(vec/norm(vec)),'descend');
vec = (vec1+vec2)/2;
vec = vec/norm(vec);
vec_arr = [vec_arr; vec kk i1 j1 nu_ang(iii,8:11);];
end

% generate stereograph colormap for visualization
color = [1 0 0;0 1 0; 0 0 1;];
[xyz,~]=spheretri(1e6);
ind = xyz(:,1)>=0 & xyz(:,2)>=0 & xyz(:,3)>=0;xyz=xyz(ind,:);
ind = xyz(:,1)>=xyz(:,2) & xyz(:,2)>=xyz(:,3);xyz=xyz(ind,:);

err = 1e-10;
[~,ind]=max(xyz(:,1));d1 = sqrt(sum((xyz-xyz(ind,:)).^2,2))*1.1;
[~,ind]=max(xyz(:,2));d2 = sqrt(sum((xyz-xyz(ind,:)).^2,2))*1.4;
[~,ind]=max(xyz(:,3));d3 = sqrt(sum((xyz-xyz(ind,:)).^2,2)); 
dis = [1./max(d1,err) 1./max(d2,err) 1./max(d3,err)];
dis = dis./max(dis,[],2);carr = dis*color;

[~,ind]=max(xyz(:,1));d1 = sqrt(sum((vec_arr(:,1:3)-xyz(ind,:)).^2,2))*1.1;
[~,ind]=max(xyz(:,2));d2 = sqrt(sum((vec_arr(:,1:3)-xyz(ind,:)).^2,2))*1.4;
[~,ind]=max(xyz(:,3));d3 = sqrt(sum((vec_arr(:,1:3)-xyz(ind,:)).^2,2)); 
dis = [1./max(d1,err) 1./max(d2,err) 1./max(d3,err)];
dis = dis./max(dis,[],2);carr1 = dis*color;

cTar = mean(xyz,1) + [0 0 4];
cPos = [0.5 4 0.5]*100;
figure(2);
subaxis(2,3,ijk1+(ijk2-1)*3,'SpacingVert',0.05,'SpacingHoriz',0.05)
scatter3(vec_arr(:,1),vec_arr(:,2),vec_arr(:,3),...
    5*ones(size(vec_arr,1),1),carr1,'filled','marker','o','linewidth',2,'sizedata',10,...
    'markeredgecolor','none');
axis image;campos(cPos + cTar);camtarget(cTar);
view(100,0);axis equal off;hold on
vec = [1 0 0;[1 1 0]/sqrt(2);[1 1 1]/sqrt(3)]';
N = 100;xx = (0:N)/N;linecal = vec(:,1)*xx+vec(:,2)*(1-xx);
for i=1:size(linecal,2)
linecal(:,i) = linecal(:,i)/norm(linecal(:,i))*1.005;
end
plot3(linecal(1,:),linecal(2,:),linecal(3,:),'k-','linewidth',2)
linecal = vec(:,1)*xx+vec(:,3)*(1-xx);
for i=1:size(linecal,2)
linecal(:,i) = linecal(:,i)/norm(linecal(:,i))*1.005;
end
plot3(linecal(1,:),linecal(2,:),linecal(3,:),'k-','linewidth',2)
linecal = vec(:,2)*xx+vec(:,3)*(1-xx);
for i=1:size(linecal,2)
linecal(:,i) = linecal(:,i)/norm(linecal(:,i))*1.005;
end
plot3(linecal(1,:),linecal(2,:),linecal(3,:),'k-','linewidth',2)
hold off
title(['\theta=' num2str(ang_cut_arr(ijk1)) '^\circ; ' num2str(dis_arr(ijk2)/10) 'nm'])
    end
end
%% plot twin boundary statistics using coincidence site lattice theory
% CSL parameters for difference boundaries
sigmaN = ...
    [  1     0 0 0 0;
       3    60 1 1 1;
       5 36.86 1 0 0;
       7 38.21 1 1 1;
       9 38.94 1 1 0;
      11 50.47 1 1 0;
     13 22.62 1 0 0;
     13 27.79 1 1 1;
      15 48.19 2 1 0;
     17 28.07 1 0 0;
     17  61.9 2 2 1;
     19 26.53 1 1 0;
     19 46.80 1 1 1;
     21 21.78 1 1 1;
     21 44.41 2 1 1;
      23 40.45 3 1 1;
     25 16.26 1 0 0;
     25 51.68 3 3 1;
     27 31.59 1 1 0;
     27 35.43 2 1 0;
     29 43.60 1 0 0;
     29 46.40 2 2 1];
indplot = [1 2 5 6 12 19 23]; % plot: sigma 1 3 9 11 19a 27a random
cryarray = importdata(['./Input/cryarray.mat']);
% plot three groups particles based on crystallinity
ind = [];
ind{1} = cryarray < 0.5;
ind{2} = cryarray > 0.5 & cryarray <0.95;
ind{3} = cryarray > 0.95;
% quantify misorientation information of selected nuclei pairs
NN = [];
for ijk = 1:3
nucleiang = cell2mat(nucleiangCell(ind{ijk})');
dis_cut = 0;
dis_range_cut = mean(valleyposition);
axis_mis_ang = 5;
nu_ang=nucleiang(abs(nucleiang(:,4)-dis_cut)<=dis_range_cut,:);
vec_arr = zeros(size(nu_ang,1),10);
for iii = 1:size(nu_ang,1)
kk=nu_ang(iii,1);i1=nu_ang(iii,2);j1=nu_ang(iii,3);
if iii == 1
data=importdata(['./Input/Nuclei_atoms_HEA_' num2str(kk) '_nanoparticle.mat']);
indGrain_combined=data.indGrain_combined;
FCCind=data.FCCindfit;
vector=data.vector;
elseif kk ~=nu_ang(iii-1,1)
data=importdata(['./Input/Nuclei_atoms_HEA_' num2str(kk) '_nanoparticle.mat']);
indGrain_combined=data.indGrain_combined;
FCCind=data.FCCindfit;
vector=data.vector;
end
vec = nu_ang(iii,8:10)*vector(:,:,i1);
vec1 = sort(abs(vec/norm(vec)),'descend');
vec = nu_ang(iii,8:10)*vector(:,:,j1);
vec2 = sort(abs(vec/norm(vec)),'descend');
vec = (vec1+vec2)/2;
vec = vec/norm(vec);
vec_arr(iii,:) = [ vec kk i1 j1 nu_ang(iii,8:11);];
end
% classify into CSL boundaries
NsigmaN = size(sigmaN,1);
Nind = zeros(NsigmaN,size(vec_arr,1));
for i=1:NsigmaN
axis_plot = sigmaN(i,3:5);
ind1 = sum(vec_arr(:,1:3).*axis_plot,2)/norm(axis_plot)>=cosd(axis_mis_ang);
ind2 = abs(vec_arr(:,10)-sigmaN(i,2)) <= min(15/sqrt(sigmaN(i,1)),1e4);
if i==1
    Nind(i,:) = ind2;
else
    Nind(i,:) = ind1 & ind2;
end
end
NN(:,ijk) = [sum(Nind,2); size(vec_arr,1) - sum(Nind,'all')];
NN(:,ijk) = NN(:,ijk)/sum(NN(:,ijk));
end
% plot CSL boundaries histogram
figure(3)
b=bar(NN(indplot,:)*100,1);
b(1).FaceColor=[0 0.5 0];
b(2).FaceColor=[0, 0.4470, 0.7410];
b(3).FaceColor=[0.6350, 0.0780, 0.1840];
box off;ylabel('Fraction (%)')
label=cell(NsigmaN,1);
for i=1:NsigmaN
    label{i}=['\Sigma' num2str(sigmaN(i,1))];
    if i~=1
    if sigmaN(i,1)==sigmaN(i-1,1)
        label{i-1}=['\Sigma' num2str(sigmaN(i,1)) 'a'];
        label{i}=['\Sigma' num2str(sigmaN(i,1)) 'b'];
    end
    end
end
label{1} = 'Merging';
label{NsigmaN+1} = 'Random';
label = label(indplot);
xticks(1:NsigmaN+1)
label = cellfun(@(x) strrep(x,'  ',','), label,'UniformOutput',false);
set(gca,'xticklabel',label);xtickangle(45);ax = axis;hold on
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1);hold off
legend(['<0.5'], ['0.5 - 0.95'],['>0.95']) % Add a legend
set(gca,'fontsize', 14,'FontName', 'Arial','fontweight','bold');legend boxoff;ax = gca;ax.BoxStyle = 'full';ax.LineWidth=1;
end