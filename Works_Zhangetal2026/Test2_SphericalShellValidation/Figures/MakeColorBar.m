function [cbar]=MakeColorBar(Position,R,G,B)
NumInput=size(Position,1);

Min=min(Position);
Max=max(Position);

NumCBar=1000;
cbar=zeros(NumCBar,3);

if NumInput>2
    Start=1;
    for i=2:NumInput
        End=((Position(i,:)-Min)/(Max-Min))*NumCBar;
        End=round(End);
        Num=End-Start+1;

        TempR=linspace(R(i-1,:),R(i,:),Num);
        TempG=linspace(G(i-1,:),G(i,:),Num);
        TempB=linspace(B(i-1,:),B(i,:),Num);
        cbar(Start:End,:)=[TempR' TempG' TempB'];
        Start=End+1;
    end
else
    TempR=linspace(R(1,:),R(2,:),NumCBar);
    TempG=linspace(G(1,:),G(2,:),NumCBar);
    TempB=linspace(B(1,:),B(2,:),NumCBar);
    cbar(1:end,:)=[TempR' TempG' TempB'];
end

cbar=cbar./255;
end