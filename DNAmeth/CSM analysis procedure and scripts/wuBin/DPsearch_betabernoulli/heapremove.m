%
% heapremove : remove the root element of the heap, and return it and the
%              new heap with the element removed.
%
% matt@lanl.gov
%
function [v,h] = heapremove(hin)
    %
    % [v,h] = heapremove(hin)
    %
    %     Remove the root of the input heap hin, and return it (v) and the
    %     new heap without it (h).
    %

    h = hin;
    v = h.tree(1);
    h.tree(1) = h.tree(h.count);
    h.count = h.count - 1;
    if (h.count == 0)
        h.tree = [];
    else
        h.tree = h.tree(1:h.count);
    end
    
    cur = 1;
    lchild = 2;
    rchild = 3;
    found = 0;
    %zhu: the following code re-arrange the order of tree after removing
    %the root.
    while (found == 0)
        numchildren = (lchild <= h.count) + (rchild <= h.count);
        
        if (numchildren == 0)
            found = 1;
        elseif (numchildren == 1)
            if (h.tree(lchild).score<h.tree(cur).score)
                tmp = h.tree(lchild);
                h.tree(lchild) = h.tree(cur);
                h.tree(cur) = tmp;
                cur = lchild;
            else
                found = 1;
            end
        else
          if h.tree(lchild).score < h.tree(rchild).score,
            if h.tree(lchild).score < h.tree(cur).score,
              tmp = h.tree(lchild);
              h.tree(lchild) = h.tree(cur);
              h.tree(cur) = tmp;
              cur = lchild;
            else
              found = 1;
            end;
          else
            if h.tree(rchild).score < h.tree(cur).score,
              tmp = h.tree(rchild);
              h.tree(rchild) = h.tree(cur);
              h.tree(cur) = tmp;
              cur = rchild;
            else
              found = 1;
            end;
          end;
        end
        
        lchild = cur*2;
        rchild = lchild+1;
    end
    
