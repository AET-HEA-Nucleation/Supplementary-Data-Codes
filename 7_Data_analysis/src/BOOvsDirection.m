clear
OPth = 0.5;
name{1}='HEA181228t1';name{2}='HEA181229t1x';name{3}='HEA181229t2';name{4}='HEA190119t1x';
name{5}='HEA190309t1x';name{6}='HEA190309t2x';name{7}='HEA190311t2x';name{8}='HEA190401t1x';
name{9}='HEA190401t2x';name{10}='HEA190408t2';name{11}='HEA190408t3';name{12}='HEA200224t1';
name{13}='HEA200302t1x';name{14}='HEA200302t2';name{15}='HEA200307t1x';name{16}='HEA200307t2x';
name{17}='HEA200920t1x';name{18}='HEA200920t2x';name{19}='HEA200924t1x';name{20}='aHEA181027t1';
name{21}='aHEA181228t2';name{22}='aHEA190119t2';name{23}='aHEA190406t1';name{24}='aHEA190406t2';
name{25}='aHEA190406t3';name{26}='aHEA190408t1';

plt=false;
list=importdata('../Analysis/2_BOP/list_arr.mat');
xx=-10:0.1:10;
order_arr123=[];
order_arr_all=[];
count=0;
for kk=1
    ii=num2str(list(kk));
outpath=['../Analysis/6_Nuclei/OPth' num2str(OPth) '/'  name{str2num(ii)} '/'];
outpath1=['../Analysis/2_BOP/' name{str2num(ii)} '/'];
% load([outpath 'NNucleiCluster.mat'])
data=importdata([outpath1 'BOPfcc0p5.mat']);
model_ori=data.model;
order=data.order;
data=importdata([outpath 'indGrain_combined_IniParaFit.mat']);
indGrain_combined=data.indGrain_combined;
vector=data.vectorfit;
data=importdata([outpath 'indGrain.mat']);
indGrain=data.indGrain;
indGrain_type = unique(indGrain(:,4));
indGrain_type(1)=[];
ngauss=4;

% N=1;
% ind=[];
% for j=1:size(NNucleiCluster{N},1)
% for i=NNucleiCluster{N}(j,:)
%     ind=[ind;find(indGrain(:,4)==i)];
% end
% end
% nucleiInd=NNucleiCluster{N}(5,1);

% nucleiInd=ind;
for nucleiInd=indGrain_type'
ind=unique(indGrain_combined(indGrain_combined(:,4)==nucleiInd,6));
if ind==-1
    continue
end
vec1=vector(:,1,ind);
vec2=vector(:,2,ind);
vec3=vector(:,3,ind);
center=indGrain(nucleiInd,1:3);
count=count+1;

for flg=1:3
if flg==1 %<100>
    NN=3;
elseif flg==2 %<110>
    NN=6;
elseif flg==3 %<111>
    NN=4;
end
order_vec_arr=zeros(NN,length(xx));

for i=1:NN
    if flg==1
        vec_arr=[vec1 vec2 vec3];
    elseif flg==2
        vec_arr=[(vec1+vec2)/sqrt(2) (vec1-vec2)/sqrt(2) (vec1+vec3)/sqrt(2) (vec1-vec3)/sqrt(2) (vec2+vec3)/sqrt(2) (vec2-vec3)/sqrt(2)];
    elseif flg==3
        vec_arr=[(vec1+vec2+vec3)/sqrt(3) (vec1-vec2+vec3)/sqrt(3) (-vec1+vec2+vec3)/sqrt(3) (-vec1-vec2+vec3)/sqrt(3)];
    end
vec=vec_arr(:,i);


xxvec=vec*xx;
pos_arr=center+xxvec';

model_nuclei=indGrain(indGrain(:,4)==nucleiInd,1:3);


shp_nuclei = alphaShape(model_nuclei(:,1),model_nuclei(:,2),model_nuclei(:,3),4);
% in_nuclei = inShape(shp_nuclei,pos_arr(:,1),pos_arr(:,2),pos_arr(:,3));

% pos_arr=pos_arr(in_nuclei,:);
% xx=xx(in_nuclei);

local_ind=sqrt(sum((model_ori-center').^2,1))<max(xx(:))+ngauss;
local_model=model_ori(:,local_ind);
local_order=order(local_ind);

% interpOBJ=rbfcreate(model_nuclei', indGrain(indGrain(:,4)==nucleiInd,5)','RBFFunction', 'gaussian', 'RBFConstant', ngauss);
interpOBJ=rbfcreate(local_model, local_order,'RBFFunction', 'linear', 'RBFConstant', ngauss);
order_vec = rbfinterp(pos_arr', interpOBJ);
order_vec_arr(i,:)=order_vec;

if plt
figure(10)
scatter3(indGrain(indGrain(:,4)==nucleiInd,1),indGrain(indGrain(:,4)==nucleiInd,2),indGrain(indGrain(:,4)==nucleiInd,3),'filled');
hold on
% scatter3(pos_arr(:,1),pos_arr(:,2),pos_arr(:,3),'k.');
plot(shp,'EdgeColor','none','FaceAlpha',0.1)
hold off
axis image off
view(35,40)
figure(11)
plot(xx,order_vec)
xlabel('distance (A)')
ylabel('BOO')
end
end
ind=find(xx==0);
xx0=xx(ind:end);
% order_arr0=mean(order_vec_arr,1);
% order_arr0=(order_arr0(ind:end)+order_arr0(ind:-1:1))/2;


order_arr0=[order_vec_arr(:,ind:end);order_vec_arr(:,ind:-1:1)];
ordercount=ones(size(order_arr0));
for j=1:NN*2
    [~,ind0]=min(order_arr0(j,:));
    if ind0<length(xx0)
        order_arr0(j,ind0:end)=0;
        ordercount(j,ind0:end)=0;
    end
end
sumordercount=sum(ordercount,1);
sumordercount(sumordercount==0)=1;
order_arr123(flg,:)=sum(order_arr0,1)./sumordercount;
order_arr_all(count,:,flg)=order_arr123(flg,:);
end
end
end
figure(1)
plot(xx0,order_arr123)
xlabel('distance (A)')
ylabel('BOO')
legend('<100>','<110>','<111>')

ind=find(xx==0);
xx0=xx(ind:end);

plotorder(1,:)=mean(order_arr_all(:,:,1),1);
plotorder(2,:)=mean(order_arr_all(:,:,2),1);
plotorder(3,:)=mean(order_arr_all(:,:,3),1);
figure(2)
plot(xx0,plotorder)
xlabel('distance (A)')
ylabel('BOO')
legend('<100>','<110>','<111>')
%% based on interface method
clear
OPth = 0.5;
name{1}='HEA181228t1';name{2}='HEA181229t1x';name{3}='HEA181229t2';name{4}='HEA190119t1x';
name{5}='HEA190309t1x';name{6}='HEA190309t2x';name{7}='HEA190311t2x';name{8}='HEA190401t1x';
name{9}='HEA190401t2x';name{10}='HEA190408t2';name{11}='HEA190408t3';name{12}='HEA200224t1';
name{13}='HEA200302t1x';name{14}='HEA200302t2';name{15}='HEA200307t1x';name{16}='HEA200307t2x';
name{17}='HEA200920t1x';name{18}='HEA200920t2x';name{19}='HEA200924t1x';name{20}='aHEA181027t1';
name{21}='aHEA181228t2';name{22}='aHEA190119t2';name{23}='aHEA190406t1';name{24}='aHEA190406t2';
name{25}='aHEA190406t3';name{26}='aHEA190408t1';

plt=false;
list=importdata('../Analysis/2_BOP/list_arr.mat');
xx=-10:0.1:10;

BOO_dis_arr=cell(3,1);
count=0;
for kk=1
ii=num2str(list(kk));
outpath=['../Analysis/6_Nuclei/OPth' num2str(OPth) '/'  name{str2num(ii)} '/'];
outpath1=['../Analysis/2_BOP/' name{str2num(ii)} '/'];

data=importdata([outpath1 'BOPfcc0p5.mat']);
model_ori=data.model;
order=data.order;
data=importdata([outpath 'indGrain_combined_IniParaFit.mat']);
indGrain_combined=data.indGrain_combined;
vector=data.vectorfit;
data=importdata([outpath 'indGrain.mat']);
indGrain=data.indGrain;
indGrain_type = unique(indGrain(:,4));
indGrain_type(1)=[];
ngauss=4;
max_radius=20;


for nucleiInd=indGrain_type'
    count=count+1;
ind=unique(indGrain_combined(indGrain_combined(:,4)==nucleiInd,6));
if ind==-1
    continue
end

modelcrystal=indGrain_combined(indGrain_combined(:,4)==nucleiInd,1:3)';
modelamorphous=model_ori(:,order<0.5);
model=[modelamorphous, modelcrystal];
order_sub=[order(order<0.5) indGrain_combined(indGrain_combined(:,4)==nucleiInd,5)'];

% modelcrystal=model_ori(:,order>0.5);
% model=model_ori;
% order_sub=order;

center=indGrain_combined(nucleiInd,1:3)';
dis=model-center;
dis=sqrt(sum(dis.^2,1));
model=model(:,dis<=max_radius);
order_sub=order_sub(dis<=max_radius);

% Method : direct calculation of atom distance to the interface
shp = alphaShape(modelcrystal(1,:)',modelcrystal(2,:)',modelcrystal(3,:)',ngauss);
facets = boundaryFacets(shp);
ind_facets=unique(facets(:));

distances=zeros(5,size(model,2));
distances(5,:)=order_sub;
for i=1:size(model,2)
    dif=model(:,i)-modelcrystal;
    dis=sqrt(sum(dif.^2,1));
    [dis_to_atom,ind_to_atom]=min(dis);
    if dis_to_atom==0 % this atom is a crystalline atom
        if inShape(shp,model(1,i),model(2,i),model(3,i)) % this atom is within alpha shape
            if sum(ind_to_atom==ind_facets) % this atom is on interface
                distances(1:3,i)=zeros(3,1);
                distances(4,i)=-1;
            else % this atom is not on interface
                % find nearest interface atom
                dif=model(:,i)-modelcrystal(:,ind_facets);
                dis=sqrt(sum(dif.^2,1));
                [dis_to_atom,ind_to_atom]=min(dis);
                temp=sum(facets==ind_facets(ind_to_atom),2)~=0;
                facets_sub=facets(temp,:);
                dis_temp=zeros(1,size(facets_sub,1));
                for j=1:length(dis_temp)
                    vec1=modelcrystal(:,facets_sub(j,2))-modelcrystal(:,facets_sub(j,1));
                    vec2=modelcrystal(:,facets_sub(j,3))-modelcrystal(:,facets_sub(j,1));
                    vec3=model(:,i)-modelcrystal(:,facets_sub(j,1));
                    nk=cross(vec1,vec2);
                    nk=nk./norm(nk);
                    height=sum(vec3.*nk);
                    vec3_in_plane=vec3-height*nk;
                    shp_2D=alphaShape([0 vec1(1) vec2(1)]',[0 vec1(2) vec2(2)]');
                    if inShape(shp_2D,vec3_in_plane(1),vec3_in_plane(2))
                        dis_temp(j)=abs(height);
                    else
                        dis_temp(j)=1e10;
                    end
                end
                if dis_temp<dis_to_atom
                    distances(1:3,i)=vec3.*nk;
                else
                    distances(1:3,i)=dif(:,ind_to_atom);
                end
                distances(4,i)=-1;
%                 distances(i)=-min([dis_temp dis_to_atom]);
            end
        else % this atom is an isolated atom
            distances(1:3,i)=zeros(3,1);
            distances(4,i)=-1;
        end
    else
                temp=sum(facets==ind_to_atom,2)~=0;
                facets_sub=facets(temp,:);
                dis_temp=zeros(1,size(facets_sub,1));
                for j=1:length(dis_temp)
                    vec1=modelcrystal(:,facets_sub(j,2))-modelcrystal(:,facets_sub(j,1));
                    vec2=modelcrystal(:,facets_sub(j,3))-modelcrystal(:,facets_sub(j,1));
                    vec3=model(:,i)-modelcrystal(:,facets_sub(j,1));
                    nk=cross(vec1,vec2);
                    nk=nk./norm(nk);
                    height=sum(vec3.*nk);
                    vec3_in_plane=vec3-height*nk;
                    shp_2D=alphaShape([0 vec1(1) vec2(1)]',[0 vec1(2) vec2(2)]');
                    if inShape(shp_2D,vec3_in_plane(1),vec3_in_plane(2))
                        dis_temp(j)=abs(height);
                    else
                        dis_temp(j)=1e10;
                    end
                end
                if dis_temp<dis_to_atom
                    distances(1:3,i)=vec3.*nk;
                else
                    distances(1:3,i)=dif(:,ind_to_atom);
                end
                distances(4,i)=1;
%                 distances(i)=min([dis_temp dis_to_atom]);
    end
end
% load vector of this nucleus
vec1=vector(:,1,ind);
vec2=vector(:,2,ind);
vec3=vector(:,3,ind);

order_dis=cell(3,1);

for flg=1:3
if flg==1 %<100>
    NN=3;
elseif flg==2 %<110>
    NN=6;
elseif flg==3 %<111>
    NN=4;
end
order_vec_arr=zeros(NN,length(xx));

for i=1:NN
    if flg==1
        vec_arr=[vec1 vec2 vec3];
    elseif flg==2
        vec_arr=[(vec1+vec2)/sqrt(2) (vec1-vec2)/sqrt(2) (vec1+vec3)/sqrt(2) (vec1-vec3)/sqrt(2) (vec2+vec3)/sqrt(2) (vec2-vec3)/sqrt(2)];
    elseif flg==3
        vec_arr=[(vec1+vec2+vec3)/sqrt(3) (vec1-vec2+vec3)/sqrt(3) (-vec1+vec2+vec3)/sqrt(3) (-vec1-vec2+vec3)/sqrt(3)];
    end
vec=vec_arr(:,i);

% dis=abs(sum(vec.*distances(1:3,:),1)).*distances(4,:);

dis=sqrt(sum(distances(1:3,:).^2,1)).*distances(4,:);

order_dis{flg}=[order_dis{flg} [dis; distances(5,:)]];
end
end

for flg=1:3
step=0.5;
xx=-5:step:10;
n=0;
nn=zeros(size(xx));
BOO_dis=zeros(size(xx));
for x=xx
    n=n+1;
    ind=abs(order_dis{flg}(1,:)-x)<step/2;
    nn(n)=sum(ind);
    BOO_dis(n)=mean(order_dis{flg}(2,ind));
end
BOO_dis_arr{flg}(count,:)=BOO_dis;

figure(11)
subplot(1,3,flg)
plot(xx,BOO_dis,'o-','linewidth',2)
ylabel('average BOO','fontsize',14)
xlabel('distance','fontsize',14)
end
% pause
end
end
%

for flg=1:3
    BOO_dis_arr{flg}(isnan(BOO_dis_arr{flg}))=0;
    count=sum(BOO_dis_arr{flg}~=0,1);
    
figure(12)
subplot(1,3,flg)
plot(xx,sum(BOO_dis_arr{flg},1)./count,'o-','linewidth',2)
ylabel('average BOO','fontsize',14)
xlabel('distance','fontsize',14)
end
