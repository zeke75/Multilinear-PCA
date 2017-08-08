function [tUs, odrIdx, TXmean, Wgt]  = MPCASYNDATA(TX,gndTX,testQ,maxK)
% MPCA for synthetic data
%
% Example: MPCASYNDATA(db1,-1,75,5)
%

%TX: (N+1)-dimensional tensor of Tensor Sample Dimension x NumSamples
N=ndims(TX)-1;%The order of samples.
IsTX=size(TX);
Is=IsTX(1:N);%The dimensions of the tensor
numSpl=IsTX(N+1);%Number of samples

%%%%%%%%%%%%%Zero-Mean%%%%%%%%%%
TXmean=mean(TX,N+1);%The mean
TX=TX-repmat(TXmean,[ones(1,N), numSpl]);%Centering

%The full projection for initialization
Qs=ones(N,1)*testQ; %%threshold
Us=cell(N,1);
tUs=cell(N,1);
Lmds=cell(N,1);

for n=1:N
    In=Is(n);Phi=zeros(In,In);
    for m=1:numSpl
        switch N
            case 2
                Xm=TX(:,:,m);
            case 3
                Xm=TX(:,:,:,m);
            case 4
                Xm=TX(:,:,:,:,m);
            otherwise
                error('Order N not supported.')
        end
        tX=tensor(Xm);
        tXn=tenmat(tX,n);
        Xn=tXn.data;
        Phi=Phi+Xn*Xn';
    end
    [Un,Lmdn]=eig(Phi);
    Lmd=diag(Lmdn);
    [stLmd,stIdx]=sort(Lmd,'descend');
    Us{n}=Un(:,stIdx); %%Un consists of eigvectors corresponding to significant eigvalues 
    tUs{n}=Us{n}';%%transpose
    Lmds{n}=Lmd(stIdx);
end

figure
plot([1:length(Lmds{1})],Lmds{1},'s-',[1:length(Lmds{2})],Lmds{2},'x-',[1:length(Lmds{3})],Lmds{3},'o-')
grid on

%Cumulative distribution of eigenvalues
cums=cell(N,1);
for n=1:N
    In=length(Lmds{n});
    cumLmds=zeros(In,1);
    Lmd=Lmds{n};
    cumLmds(1)=Lmd(1);
    for in=2:In
        cumLmds(in)=cumLmds(in-1)+Lmd(in);
    end
    cumLmds=cumLmds./sum(Lmd);
    cums{n}=cumLmds;
end

figure
plot([1:length(cums{1})],cums{1},'s-',[1:length(cums{2})],cums{2},'x-',[1:length(cums{3})],cums{3},'o-')
grid on

%MPCA Iterations
IPM=cell(3,1); %% 3 choices of initialization

tPs=cell(N,1);
pUs=cell(N,1);    
%%%%%%%%%%%%%Determine Rn, the dimension of projected space%%%%
for n=1:N
    cum=cums{n};
    idxs=find(cum>=Qs(n)/100);%%Q-based method
    Ps(n)=idxs(1);%%dimension kept
    tUn=tUs{n};
    tPn=tUn(1:Ps(n),:);
    tPs{n}=tPn;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IPM{1}=tPs; %full projection truncation

for n=1:N
    tIm{n} = eye(Ps(n),Is(n));
end
IPM{2}=tIm; %truncated identity matrices

for n=1:N
    rPm{n} = 10*rand(Ps(n),Is(n));
end
IPM{3}=rPm; %random matrices

PsiY=cell(3,1);
for init=1:3
    for iK=1:maxK
        for n=1:N
            In=Is(n);%%In is full projection dimension
            Phi=double(zeros(In,In));
            for m=1:numSpl
                switch N
                    case 2
                        Xm=TX(:,:,m);
                    case 3
                        Xm=TX(:,:,:,m);
                    case 4
                        Xm=TX(:,:,:,:,m);
                    otherwise
                        error('Order N not supported.')
                end
                tX=ttm(tensor(Xm),IPM{init},-n); %%n-way multiply note "-n" multiply except n-mode
                tXn=tenmat(tX,n);
                Xn=tXn.data;
                Phi=Phi+Xn*Xn';
            end
            Pn=Ps(n);
            Phi=double(Phi);
            if Pn<In
                option=struct('disp',0);
                [pUs{n},pLmdn]=eigs(Phi,Pn,'lm',option);
                pLmds{n}=diag(pLmdn);
            else
                [pUn,pLmdn]=eig(Phi);
                pLmd=diag(pLmdn);
                [stLmd,stIdx]=sort(pLmd,'descend');
                pUs{n}=pUn(:,stIdx(1:Pn));
                pLmds{n}=pLmd(stIdx(1:Pn));
            end
            %tPs{n}=pUs{n}';
            IPM{init}{n}=pUs{n}';

        end

        Y=ttm(tensor(TX),IPM{init},1:N);
        Y=Y.data;
        Psi=0;

        for m=1:numSpl    
            switch N
                case 2
                    Ym=Y(:,:,m);
                case 3
                    Ym=Y(:,:,:,m);
                case 4
                    Ym=Y(:,:,:,:,m);
                otherwise
                    error('Order N not supported.')
            end    
            Psi=Psi+norm(tensor(Ym))^2;
        end
        PsiY{init}=[PsiY{init},Psi]     
    end
end
Us=pUs;
tUs=tPs;
Lmds=pLmds;
Is=Ps;

figure
plot([1:maxK],PsiY{1},'s-',[1:maxK],PsiY{2},'x-',[1:maxK],PsiY{3},'o-')
grid on


