function u=getMC1(kw,bdot,mmax)

    u=-kw*bdot;

    for j=1:length(u)

        if abs(u(j))>mmax

        u(j)=sign(u(j))*mmax;

    end

end