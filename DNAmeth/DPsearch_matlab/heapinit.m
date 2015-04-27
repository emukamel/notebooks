%
% heapinit
%
% matt@lanl.gov
% 
function h = heapinit()
    %
    % h = heapinit
    %
    %     Return a heap h that is empty.  This must be called before the
    %     heapinsert and heapremove routines.
    %
    h.count = 0;
    h.tree = [];
    
