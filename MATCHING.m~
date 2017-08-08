function idenrate=MATCHING(Prbset,rank,Hz,k)

% Matching algorithm for recognition
% The rtesult is not good.


load Gal

switch Prbset
    case 1
        load PrbA
    case 2
        load PrbB
    case 3
        load PrbC
    case 4
        load PrbD
    case 5
        load PrbE
    case 6
        load PrbF
    case 7
        load PrbG
end

[tUs,odrIdx,TXmean,Wgt,vecYps]  = MPCA(fea3D,gnd,97,1);
[tUsg,odrIdxg,TXmeang,Wgtg,vecYpsg]  = MPCA(fea3Dg,gndg,97,1);

odrIdx=odrIdx(1:Hz);%Take only the first a few
vecYps=vecYps(odrIdx,:);%Input to LDA

odrIdxg=odrIdxg(1:Hz);%Take only the first a few
vecYpsg=vecYpsg(odrIdxg,:);%Input to LDA


classLabel = unique(gnd);
nClass = length(classLabel);%Number of classes
ClsIdxs=cell(nClass);
Ns=zeros(nClass,1);
for i=1:nClass
    ClsIdxs{i}=find(gnd==classLabel(i));
    Ns(i)=length(ClsIdxs{i}); %Number of samples in ith class
end

classLabelg = unique(gndg);
ClsIdxsg=cell(71);
Nsg=zeros(71,1);
for i=1:71
    ClsIdxsg{i}=find(gndg==classLabelg(i));
    Nsg(i)=length(ClsIdxsg{i}); %Number of samples in ith class
end

%Compute score between class p in probe set and class g in gallery set
s=zeros(nClass,71); %score s(p,g)
for p=1:nClass %probe set
    for g=1:71 %gallery set
        fv=vecYps(:,ClsIdxs{p});
        fvg=vecYpsg(:,ClsIdxsg{g});
        d=zeros(Ns(p),Nsg(g));
        for i=1:Ns(p)
            for j=Nsg(g)
                d(i,j)=distance(fv(:,i),fvg(:,j),1,k);
            end
        end
        s(p,g)=mean(-min(d'))+mean(-min(d));
    end
end


%Matching
n=0;
for p=1:nClass
    [a,order]=sort(s(p,:),'descend');
    for i=1:rank
        if order(i)==p
            n=n+1;
        end
    end    
end
idenrate=n/nClass;
        
        
        
    

        
        
        
        
            
 
