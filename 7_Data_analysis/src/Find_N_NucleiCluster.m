clear
for kk=50:53
    kk
% list=importdata('../Analysis/2_BOP/list_arr.mat');
list=1:53;
ii=num2str(list(kk));
OPth = 0.5;
% ii=num2str(kk);

N=10; % find up to 10 nuclei clusters

name{1}='HEA181228t1';name{2}='HEA181229t1x';name{3}='HEA181229t2';name{4}='HEA190119t1x';
name{5}='HEA190309t1x';name{6}='HEA190309t2x';name{7}='HEA190311t2x';name{8}='HEA190401t1x';
name{9}='HEA190401t2x';name{10}='HEA190408t2';name{11}='HEA190408t3';name{12}='HEA200224t1';
name{13}='HEA200302t1x';name{14}='HEA200302t2';name{15}='HEA200307t1x';name{16}='HEA200307t2x';
name{17}='HEA200920t1x';name{18}='HEA200920t2x';name{19}='HEA200924t1x';name{20}='aHEA181027t1';
name{21}='aHEA181228t2';name{22}='aHEA190119t2';name{23}='aHEA190406t1';name{24}='aHEA190406t2';
name{25}='aHEA190406t3';name{26}='aHEA190408t1';
name{50}='PdAt1';name{51}='PdA350t1';name{52}='PdAt2';name{53}='PdA350t2';

peakarray=importdata('../Analysis/1_RDF/peakposition.mat');
valleyposition=importdata('../Analysis/1_RDF/valleyposition.mat');

inpath=['../Analysis/2_BOP/' name{str2num(ii)} '/'];
% outpath=['../Analysis/6_Nuclei/' name{str2num(ii)} '/'];
outpath=['../Analysis/6_Nuclei/OPth' num2str(OPth) '/'  name{str2num(ii)} '/'];

data=importdata([inpath 'BOPfcc0p5.mat']);
scaled_SROP=data.order;
model=double(data.model);
atoms=data.atom;

data=importdata([outpath 'indGrain.mat']);
indGrain=data.indGrain;
indGrain_type = unique(indGrain(:,4));
indGrain_type(1)=[];

NNucleiCluster=cell(N,1);
%% find nuclei neighbor
neighbor=[];
count=length(indGrain_type);
for i=1:count-1
    i
    for j=i+1:count
        ind1=find(indGrain(:,4)==indGrain_type(i));
        ind2=find(indGrain(:,4)==indGrain_type(j));
        model1=model(:,ind1);
        model2=model(:,ind2);
        
        flg1=false;
        for k=1:length(ind1)
            dif=model2-repmat(model1(:,k),[1 length(ind2)]);
            dis=sqrt(sum(dif.^2,1));
            if sum(dis<=valleyposition(str2num(ii)))
                flg1=true;
                break
            end
        end
        if flg1
            neighbor=[neighbor;indGrain_type(i) indGrain_type(j) i j;];
        end
    end
end
%% find N nuclei clusters
if isempty(neighbor)
    NNucleiCluster{1}=indGrain_type(:);
else
counted=[];
for i=1:count
    i
    if sum(i==counted)~=0
        continue
    end
    ind=find(sum(neighbor(:,3:4)==i,2)~=0);
    temp=neighbor(ind,3:4);
    temp=unique(temp(:));
    ind_cluster=temp;
    j=2;
    while j<=length(ind_cluster)
        ind=find(sum(neighbor(:,3:4)==ind_cluster(j),2)~=0);
        temp=neighbor(ind,3:4);
        temp=unique(temp(:));
        for k=temp'
            if sum(k==ind_cluster)==0
                ind_cluster=[ind_cluster;k];
            end
        end
        j=j+1;
    end
    Nsize=length(ind_cluster);
    if Nsize==0
        NNucleiCluster{1}=[NNucleiCluster{1};indGrain_type(i)];
    elseif Nsize<=N
        NNucleiCluster{Nsize}=[NNucleiCluster{Nsize};indGrain_type(ind_cluster)'];
    end
    counted=[counted;ind_cluster];
end
end
save([outpath 'NNucleiCluster.mat'],'NNucleiCluster','neighbor')
end
%%
N=4;
figure(11);clf;hold on;
for j=1:size(NNucleiCluster{N},1)
for i=NNucleiCluster{N}(j,:)
    figure(11)
    scatter3(indGrain(indGrain(:,4)==i,1),indGrain(indGrain(:,4)==i,2),indGrain(indGrain(:,4)==i,3),'filled','markeredgecolor','k');
    axis image
    title(['Total atoms: ' num2str(sum(indGrain(:,4)==i))])
%     pause
end
end
view(35,40)
