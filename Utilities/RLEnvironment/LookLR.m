function [val,valid,w] = LookLR(s,w,row,col,actName)
global world2D;
global tempworld2D;

valid = 0;
w.limbTeleLR=0;

if strcmp(s.wrld.resetType,'BottomTop_InfLR')
    valid = 1;
    switch actName
        case {'LEFT','LEFTBUTTON'} % move left
            if col-1 >= 2 & col-1 <= size(world2D,2)-1
                val = tempworld2D(row,col-1);
            elseif col-1 < 2
                val = tempworld2D(row,end-1);
                w.limbTeleLR=1;
            end
        case {'STAY','BUTTON'} % stay still
            if col >= 2 & col <= size(world2D,2)-1
                val = tempworld2D(row,col);
            end
        case {'RIGHT','RIGHTBUTTON'} % move right
            if col+1 >= 2 & col+1 <= size(world2D,2)-1
                val = tempworld2D(row,col+1);
            elseif col+1 > size(world2D,2)-1
                val = tempworld2D(row,2);
                w.limbTeleLR=-1;
            end
    end
    
else
    
    switch actName
        case 'LEFT' % move left
            if col-1 >= 1 & col-1 <= size(world2D,2)
                val = tempworld2D(row,col-1);
                valid = 1;
            end
        case {'STAY','BUTTON'} % stay still
            if col >= 1 & col <= size(world2D,2)
                val = tempworld2D(row,col);
                valid = 1;
            end
        case 'RIGHT' % move right
            if col+1 >= 1 & col+1 <= size(world2D,2)
                val = tempworld2D(row,col+1);
                valid = 1;
            end
    end
    
end


end