function Syndata=build_dataset(f)

% Construct the synthetic data (Section V.A in the paper)
%
% f controls the eigenvalue distribution

Dim=[30,20,10];
c=cell(3,1);
for i=1:3
    [U,S,V]=svd(normrnd(0,1,[Dim(i),Dim(i)]));
    c{i}=U; %%Projection matrix C
end

Syndata=[];
for m=1:100
    for i1=1:30
        for i2=1:20
            for i3=1:10
                 B(i1,i2,i3)=normrnd(0,1)*((6000/(i1*i2*i3))^f); %%core tensor B
            end
        end
    end   
    D=normrnd(0,0.01,Dim); %%noise tensor D
    A=ttm(tensor(B),c)+D;
    Syndata(:,:,:,m)=A;
end
    