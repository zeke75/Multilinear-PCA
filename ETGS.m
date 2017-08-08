function eigentensor=ETGS(TX,gndTX,testQ,maxK,k)

% Obtain the kth discrimitive eigentensor
%
% Example: %ETGS(fea3D,gnd,97,1,1), %ETGS(fea3D,gnd,97,1,2),...
%

[tUs,odrIdx,TXmean,Wgt] = MPCA(TX,gndTX,testQ,maxK);

sU1=size(tUs{1});
Is(1)=sU1(1);
sU2=size(tUs{2});
Is(2)=sU2(1);
sU3=size(tUs{3});
Is(3)=sU3(1);

%Find the position of feature vector
p=zeros(3,1);
%
p(1)=mod(odrIdx(k),Is(1));
if p(1)==0 
    p(1)=Is(1);
end

p(2)=ceil(odrIdx(k)/(Is(1)*Is(3)));

p(3)=mod(ceil(odrIdx(k)/(Is(1))),Is(3));
if p(3)==0 
    p(3)=Is(3);
end

U1=tUs{1};
u1=U1(p(1),:);
U2=tUs{2};
u2=U2(p(2),:);
U3=tUs{3};
u3=U3(p(3),:);

etg=ktensor({u1',u2',u3'});
etg=tensor(etg);
etg=tenmat(etg,1);
etg=etg.data
ma=max(max(etg));
mi=min(min(etg));
eigentensor=(etg-mi)/(ma-mi);

imshow(eigentensor)
