function [row,column] = K_index(IEN_all,Grid_size) 
keyCount = 0;
row = [];
column = [];
%Kindex = [];
%BCindex = [];
%BCeleIndex = [];
for ele = 1:Grid_size 
    IENall = IEN_all(ele,:);
    nssl = length(IENall);
    for i = 1:nssl
        rowTemp = IENall(i);
         for j = 1:nssl           
            columnTemp = IENall(j);
            keyCount = keyCount+1; %its a new key
            row = [row, rowTemp];
            column = [column, columnTemp];
            
%             if ismember(ele,BCB)
%                 BCeleIndex = [BCeleIndex,ele];
%                 BCindex = [BCindex,keyCount];
%             end
            
         end
    end
end

% rowIndex = zeros(length (IBC));
% columnIndex = zeros(length (IBC));
% for i = 1: length (IBC)
%     rowI = find(row==IBC(i));
%     columnI = find(column==IBC(i));
%     for j = 1:length(rowI)
%         rowIndex(i,j) = rowI(j);
%     end
%     for j = 1:length(columnI)
%         columnIndex(i,j) = columnI(j);
%     end
% end
% 
% BCB_Kindex = [];
% for i = 1:length(BCB)
% %for i = 1:1   
%     bcele = BCB(i);
%     indexBCBindex = bcele == BCeleIndex;
%     BCBindex = BCindex(indexBCBindex);
%     BCB_Kindex = [BCB_Kindex,BCBindex];
% end

end