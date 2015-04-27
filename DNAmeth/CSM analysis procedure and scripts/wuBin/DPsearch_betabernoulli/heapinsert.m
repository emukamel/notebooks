%
% heapinsert
%
% matt@lanl.gov
%
function hout = heapinsert(hin,val)
    %
    % hout = heapinsert(hin,val)
    %
    %    return the heap created by inserting the value val into the
    %    existing heap hin.
    %
    if (hin.count == 0)
        hout.count = 1;
        hout.tree = [val];
    else
        hout.count = hin.count+1;
        hout.tree = [hin.tree val]; % a matrix of structures. 
    end
    
    cur = hout.count; 
    parent = floor(cur/2);
    found = 0;
    
    while (found == 0) % zhu: the following chunk of code is to make sure the order of the tree are by increasing order of scores.
        if (parent == 0)
            found = 1;
        else
            % zhunote: move the order of the tree so that the smaller scores goes to the left.  
            if (hout.tree(parent).score > hout.tree(cur).score)
                tmp = hout.tree(parent);
                hout.tree(parent) = hout.tree(cur);
                hout.tree(cur) = tmp;
                cur = parent;
            else
                found = 1;
            end
        end
        
        parent = floor(cur/2);
    end
    
    
