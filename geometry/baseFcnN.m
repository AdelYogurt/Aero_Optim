function N=baseFcnN(u_list,u_x,i,k)
% base function of BSpline curve
%
if k == 0
    if ((u_list(i) <= u_x) && (u_x <= u_list(i+1)))
        if any(u_list == u_x) && u_x ~= u_list(1) && u_x ~= u_list(end)
            N=0.5;
        else
            N=1;
        end
    else
        N=0;
    end
else
    if u_list(i+k) == u_list(i)
        A=0;
    else
        A=(u_x-u_list(i))/(u_list(i+k)-u_list(i));
    end

    if u_list(i+k+1) == u_list(i+1)
        B=0;
    else
        B=(u_list(i+k+1)-u_x)/(u_list(i+k+1)-u_list(i+1));
    end

    N=A*baseFcnN(u_list,u_x,i,k-1)+B*baseFcnN(u_list,u_x,i+1,k-1);
end
end
