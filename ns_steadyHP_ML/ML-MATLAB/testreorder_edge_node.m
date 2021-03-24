function testreorder_edge_node()

Lev =[

     2     2
     5     2
     2     2
     2     2];


eV =[  1128         759           0           0           0                   
        1131        1127       10236       10237       10238                   
        1128        1131           0           0           0                
         759        1127           0           0           0];


nEv =17;


vID =[1128
        1131
         759
        1127
     0
     0
     0
     0
     0];


nR =4;

vID = updateExisting(Lev,eV,nEv,vID,nR)

end

function [vID] = updateExisting(lEv_edge,eV,nEv,vID,level)
for e = 1:4
    if nEv>2
        nEv
        lEv = lEv_edge(e,1)
        levlevel = log2(lEv-1);
        if nEv <=lEv %% all edge node are existing
            eV(e,:)
            Reorder_eV = reorder_edge_node(eV(e,:),levlevel,nEv) %% order edge node to levels
            evstart = 5+(nEv-2)*(e-1)
            vID(evstart:evstart+nEv-3) = convert_back(Reorder_eV(1:nEv-2),level,nEv) %% convert back
        else %% some or all edge node will need to be created
            if lEv>2
                
                Reorder_eV = reorder_edge_node(eV(e,:),levlevel,lEv) %% order edge node to levels
                evstart = 5+(nEv-2)*(e-1)
                temp = Reorder_eV(1:lEv-2);
                tempIndex = reorder_index(level)-3
                for index = 1:lEv-2
                    evstart+tempIndex(index)
                   vID(evstart+tempIndex(index)) = temp(index);
                end
            end
        end
    end
end
%vID
end

function Reorder_eV = reorder_edge_node(eV,level,lEv)
level
lEv
Reorder_index = reorder_index(level)
Reorder_eV = zeros(lEv-2,1);
for i = 1:lEv-2
    index = Reorder_index(i);
    Reorder_eV(i) = eV(index);
end
%Reorder_index
%Reorder_eV = Reorder_eV(3:end);
end

function converted_eV = convert_back(eV,level,lEv)
Reorder_index = reorder_index(level);
converted_eV = zeros(lEv-2,1);
for i = 1:lEv-2
    index = Reorder_index(i);
    converted_eV(index-2) = eV(i);
end
%converted_eV = converted_eV(3:end)
end

function  Reorder_index = reorder_index(level)
%% reorder edge node based on levels

temp = zeros(level,1);
for i = 2:level+1
    temp(i-1) = 2^(i-2);
end

s = 3;
start = zeros(level,1);
interval = zeros(level,1);
for i = 1:level
    s = s+temp(i);
    start(i+1) = s;
    interval(i) = 2^i; 
end
start(1) = 3;

if level>0
    count = 1;
    for i = level:-1:1
        j = 2^(i-1)+2;
        k = 2^i+1;
        
        ss = start(count);
        in = interval(count);
        
        c = 0;
        for ll = j:k
            Reorder_index(ll) = ss+in*c;
            c = c+1;
        end
        count = count+1;
    end
end

Reorder_index
Reorder_index = Reorder_index(3:end);
end

