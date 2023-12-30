
function L=baseFcnL(u_list,u_new_list,i,j,k)
% base function of increase Degree
%
if k == 0
    if ((u_list(i) <= u_new_list(j)) && (u_new_list(j) <= u_list(i+1)))
        if any(u_list == u_new_list(j)) && u_new_list(j) ~= u_list(1) && u_new_list(j) ~= u_list(end)
            L=0.5;
        else
            L=1;
        end
    else
        L=0;
    end
else
    if u_list(i+k) == u_list(i)
        A=0;
    else
        A=(u_new_list(j+k+1)-u_list(i))/(u_list(i+k)-u_list(i));
    end

    if u_list(i+k+1) == u_list(i+1)
        B=0;
    else
        B=(u_list(i+k+1)-u_new_list(j+k+1))/(u_list(i+k+1)-u_list(i+1));
    end

    L=A*baseFcnL(u_list,u_new_list,i,j,k-1)+B*baseFcnL(u_list,u_new_list,i+1,j,k-1)+baseFcna(u_list,u_new_list,i,j,k);
end
end

function a=baseFcna(u_list,u_new_list,i,j,k)
% discrete base function of BSpline curve
%
if k == 0
    if ((u_list(i) <= u_new_list(j)) && (u_new_list(j) < u_list(i+1)))
        if any(u_list == u_new_list(j)) && u_new_list(j) ~= u_list(1) && u_new_list(j) ~= u_list(end)
            a=0.5;
        else
            a=1;
        end
    else
        a=0;
    end
else
    if u_list(i+k) == u_list(i)
        A=0;
    else
        A=(u_new_list(j+k)-u_list(i))/(u_list(i+k)-u_list(i));
    end

    if u_list(i+k+1) == u_list(i+1)
        B=0;
    else
        B=(u_list(i+k+1)-u_new_list(j+k))/(u_list(i+k+1)-u_list(i+1));
    end

    a=A*baseFcna(u_list,u_new_list,i,j,k-1)+B*baseFcna(u_list,u_new_list,i+1,j,k-1);
end
end
