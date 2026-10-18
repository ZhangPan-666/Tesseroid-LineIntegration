function [CData]=GetCData(Data,ColorBar,Range)
if ~isempty(Range)
    MaxRange=max(Range);
    MinRange=min(Range);
else
    MaxRange=max(Data);MaxRange=MaxRange(1);
    MinRange=min(Data);MinRange=MinRange(1);
end

Data(Data>MaxRange)=MaxRange;
Data(Data<MinRange)=MinRange;

Data=(Data-MinRange)/(MaxRange-MinRange);

[Rows,Cols]=size(Data);
NumColor=size(ColorBar,1);

CData=zeros(Rows,Cols,3);

for i=1:Rows
    for j=1:Cols
        Position=round(Data(i,j)*NumColor);
        Position(Position>NumColor)=NumColor;
        Position(Position<1)=1;

        CData(i,j,1)=ColorBar(Position,1);
        CData(i,j,2)=ColorBar(Position,2);
        CData(i,j,3)=ColorBar(Position,3);
    end
end

end