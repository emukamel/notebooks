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
        hout.tree = [hin.tree val];
    end
    
    cur = hout.count;
    parent = floor(cur/2);
    found = 0;
    
    while (found == 0)
        if (parent == 0)
            found = 1;
        else
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
    
    
