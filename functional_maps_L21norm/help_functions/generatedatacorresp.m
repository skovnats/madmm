function [MSER]=generatedatacorresp(F0,n)
%%
FX=F0{1};
FY=F0{2};

nf=size(FX,2);
MSERX=FX;
MSERY=FY;
%% 0 - sum / 1 - multiplication
for i=1:n
    %% generate a random number of functions
    r=randint(nf,1,1);
    %% sample r random funtions
    %    ind=randint(nf,1,r);
    ind=randperm(nf);
    ind=ind(1:r);
    %%
    for j=1:(length(ind)-1)
        icur=ind(j);
        inext=ind(j+1);
        %% what to do with functions
        a=randint(2,1,1)-1;
        fX=zeros(size(FX(:,1)));
        fY=zeros(size(FY(:,1)));
        switch a
            case 0
                fX=fX+FX(:,icur)+FX(:,inext);
                fY=fY+FY(:,icur)+FY(:,inext);
            case 1
                fX=fX+FX(:,icur).*FX(:,inext);
                fY=fY+FY(:,icur).*FY(:,inext);
        end
        %%
        fX=double(logical(fX));
        fY=double(logical(fY));
    end
    %%
    fX=double(logical(fX));
    fY=double(logical(fY));
    
    %%
    MSERX=[MSERX,fX];
    MSERY=[MSERY,fY];
end
MSER{1}=MSERX;
MSER{2}=MSERY;
end