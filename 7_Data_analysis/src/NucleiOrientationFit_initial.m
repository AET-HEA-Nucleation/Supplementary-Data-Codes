function [uFCCind,vec,blenFit,ori]=NucleiOrientationFit_initial(indmodel,indmodel1,ori,vector1,blen,del)
plt=0;
% del=0.15;
vec1=vector1(:,1);
vec2=vector1(:,2);
vec3=vector1(:,3);
vec0=vector1;
angle=getAngle(vec0)'; %ZYX angle in rad

uCorner=[];
for i=1:size(indmodel1,2)
    dis=sum(abs(indmodel-indmodel1(:,i)),1);
    [dis,ind]=min(dis);
    if dis<1e-3
        uCorner=[uCorner ind];
    end
end

dif=indmodel(:,uCorner)-repmat(ori,[1 length(uCorner)]);
        
xx=round(sum(dif.*vec1,1)/blen*2)/2;
yy=round(sum(dif.*vec2,1)/blen*2)/2;
zz=round(sum(dif.*vec3,1)/blen*2)/2;

% dx=xx-sum(dif.*vec1,1)/blen;dy=yy-sum(dif.*vec2,1)/blen;dz=zz-sum(dif.*vec3,1)/blen;
% 
% dis=sqrt(dx.^2+dy.^2+dz.^2)*blen;
% uCorner=find(dis<=blen*del);

xdata.uCorner=uCorner;
xdata.xx=xx;%(uCorner);
xdata.yy=yy;%(uCorner);
xdata.zz=zz;%(uCorner);
xdata.indmodel=indmodel;
x0(1:3)=ori;x0(4:6)=angle;x0(7)=blen;
lb=x0;ub=x0;
lb(1:3)=ori-blen*del;lb(4:6)=angle-10/180*pi;lb(7)=blen-1;
ub(1:3)=ori+blen*del;ub(4:6)=angle+10/180*pi;ub(7)=blen+1;
        
opt = optimset('TolFun',1e-12);
%%%%%% use initial values and no fit with flg=0
% uFCCind=uCorner;
% para0=x0;
%%%%%%%%%%%%%%%%%
flg=1; %turn off fitting
niter=0;
while flg
    niter=niter+1;
if length(xdata.uCorner)<3
    para0=x0;
    uFCCind=xdata.uCorner;
    break
end
[para0, ~,~] = lsqcurvefit(@FCClatticeErr, x0, xdata, zeros(3,length(xdata.uCorner)), lb, ub, opt);
        
% para0=x0;

orif=para0(1:3)';
anglef=para0(4:6);
blenf=para0(7);
vecfit=getMatrix(anglef);
vecf1=vecfit(:,1);vecf2=vecfit(:,2);vecf3=vecfit(:,3);

ref=repmat(orif,[1 length(xdata.uCorner)])+(vecf1*xdata.xx+vecf2*xdata.yy+vecf3*xdata.zz)*blenf;

if plt
figure(1)
scatter3(indmodel(1,:),indmodel(2,:),indmodel(3,:),'o');hold on
scatter3(indmodel(1,xdata.uCorner),indmodel(2,xdata.uCorner),indmodel(3,xdata.uCorner),'o','filled')
scatter3(ref(1,:),ref(2,:),ref(3,:),'k.');

for i=1:length(xdata.uCorner)-1
    for j=i+1:length(xdata.uCorner)
        dis=ref(:,i)-ref(:,j);
        dis=sqrt(sum(dis.^2,1));
        if dis<3.5
            hold on
            line(ref(1,[i j]),ref(2,[i j]),ref(3,[i j]),'color','k')
        end
    end
end


hold off
axis image off
end

vec1=vecf1;vec2=vecf2;vec3=vecf3;
dref=[vec1,-vec1,vec2,-vec2,vec3,-vec3,(vec1+vec2)/2,(vec1-vec2)/2,(-vec1+vec2)/2,(-vec1-vec2)/2,...
    (-vec2+vec3)/2,(vec2-vec3)/2,(-vec2-vec3)/2,(vec2+vec3)/2,...
	(vec1+vec3)/2,(vec1-vec3)/2,(-vec1+vec3)/2,(-vec1-vec3)/2]*blenf; %for FCC

%%
uFCCind=xdata.uCorner(1);
ref_FCC=orif;
n=1;
while n<=length(uFCCind)
    ref_temp=ref_FCC(:,n)+dref;
    for i=1:18
        dif=indmodel-ref_temp(:,i);
        dis=sqrt(sum(dif.^2,1));
        ind=find(dis<=del*blenf);
        if ~isempty(ind) && ~ismember(ind,uFCCind)
            uFCCind(end+1)=ind;
            ref_FCC(:,end+1)=ref_temp(:,i);
        end
    end
    n=n+1;
end
% if the starting atom is too far from lattice, remove it
dif=indmodel(:,xdata.uCorner(1))-orif;
dis=sqrt(sum(dif.^2));
if dis>del*blenf
    uFCCind(1)=[];
    ref_FCC(:,1)=[];
end

if length(uFCCind)<length(xdata.uCorner) && niter>1 % fitted model has less atoms than that before fitting
    uFCCind=xdata.uCorner;
    para0=x0;
    break
end

if plt
figure(1)
hold on
scatter3(indmodel(1,uFCCind),indmodel(2,uFCCind),indmodel(3,uFCCind),'k*');hold off
end
if niter>10 || length(xdata.uCorner)==length(uFCCind)
    flg=0;
end
%%
ori=orif;
dif=indmodel(:,uFCCind)-repmat(ori,[1 length(uFCCind)]);
        
xx=round(sum(dif.*vec1,1)/blenf*2)/2;
yy=round(sum(dif.*vec2,1)/blenf*2)/2;
zz=round(sum(dif.*vec3,1)/blenf*2)/2;

xdata.xx=xx;
xdata.yy=yy;
xdata.zz=zz;
xdata.indmodel=indmodel;
xdata.uCorner=uFCCind;

x0=para0;
lb=x0;ub=x0;
lb(1:3)=ori-blenf*del;lb(4:6)=para0(4:6)-10/180*pi;lb(7)=blenf-1;
ub(1:3)=ori+blenf*del;ub(4:6)=para0(4:6)+10/180*pi;ub(7)=blenf+1;
if plt
pause
end
end
blenFit=para0(7);
vec=getMatrix(para0(4:6));
ori=para0(1:3)';