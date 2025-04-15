function dif = FCClatticeErr(para, xdata)
ori=para(1:3)';
angle=para(4:6);
blen=para(7);

vec=getMatrix(angle);
vec1=vec(:,1);vec2=vec(:,2);vec3=vec(:,3);

xx=xdata.xx;yy=xdata.yy;zz=xdata.zz;
indmodel=xdata.indmodel;
uCorner=xdata.uCorner;

ref=repmat(ori,[1 length(uCorner)])+(vec1*xx+vec2*yy+vec3*zz)*blen;

dif=indmodel(:,uCorner)-ref;