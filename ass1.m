function [] = ass1(num)
    if (num==1)    
        lambda = 0.8; 
        tinit = 0; 
        tend = 30;
        gp = 41;
        delx=2/(gp-1);
        x1 = -1;
        %uinitial = zeros(1:gp);
        for x=1:(gp)
            r=x1+(x-1)*delx;
            uinitial(x)=ux1(r);
        end
        delt = delx*lambda;
        it = round((tend-tinit)/delt);
        u = uinitial;
        fu = testcase1(u);
        n = length(uinitial);
        for i=1:it
            nextu(1) = 1/2*(u(n)+u(2)) + lambda/2*(fu(n) - fu(2));
            for j=2:n-1
                nextu(j) = 1/2*(u(j-1)+u(j+1)) + lambda/2*(fu(j-1) - fu(j+1));
            end
            nextu(n) = 1/2*(u(n-1)+u(1)) + lambda/2*(fu(n-1) - fu(1));
            u = nextu;
            fu = testcase1(u);
        end
        final = u;
        plot(0:40, final, 'ro', 0:40, uinitial, 'b-');
        
    else if (num == 2)
            lambda = 0.8; 
            tinit = 0; 
            tend = 4;
            gp = 41;
            delx = 2/(gp);
            x1 = -1;
            for x=1:gp
                r=x1+(x-1)*delx;
                uinitial(x)=ux2(r);
            end
            delt = delx*lambda;
            it = round((tend-tinit)/delt);
            u = uinitial;
            fu = testcase1(u);
            n = length(uinitial);
            for i=1:it
                nextu(1) = 1/2*(u(n)+u(2)) + lambda/2*(fu(n) - fu(2));
                for j=2:n-1
                    nextu(j) = 1/2*(u(j-1)+u(j+1)) + lambda/2*(fu(j-1) - fu(j+1));
                end
                nextu(n) = 1/2*(u(n-1)+u(1)) + lambda/2*(fu(n-1) - fu(1));
                u = nextu;
                fu = testcase1(u);
            end
            final = u;
            plot(0:40, final, 'ro', 0:40, uinitial, 'b-');
        else if (num == 3)
                lambda = 0.8; 
                tinit = 0; 
                tend1 = 4;
                tend2 = 40;
                gp = 601;
                delx = 2/(gp);
                x1 = -1;
                for x=1:gp
                    r=x1+(x-1)*delx;
                    uinitial(x)=ux2(r);
                end
                delt = delx*lambda;
                iux1 = round((tend1-tinit)/delt);
                iux2 = round((tend2-tinit)/delt);
                delt = delx*lambda;
                u1 = uinitial;
                fu1 = testcase1(u1);
                n = length(uinitial);
                for i=2:iux1
                    nextu1(1) = 1/2*(u1(n)+u1(2)) + lambda/2*(fu1(n) - fu1(2));
                    for j=2:n-1
                        nextu1(j) = 1/2*(u1(j-1)+u1(j+1)) + lambda/2*(fu1(j-1) - fu1(j+1));
                    end
                    nextu1(n) = 1/2*(u1(n-1)+u1(1)) + lambda/2*(fu1(n-1) - fu1(1));
                    u1 = nextu1;
                    fu1 = testcase1(u1);
                end
                final1 = u1;
                delt = delx*lambda;
                u2 = uinitial;
                fu2 = testcase1(u2);
                n = length(uinitial);
                for i=2:iux2
                    nextu2(1) = 1/2*(u2(n)+u2(2)) + lambda/2*(fu2(n) - fu2(2));
                    for j=2:n-1
                        nextu2(j) = 1/2*(u2(j-1)+u2(j+1)) + lambda/2*(fu2(j-1) - fu2(j+1));
                    end
                    nextu2(n) = 1/2*(u2(n-1)+u2(1)) + lambda/2*(fu2(n-1) - fu2(1));
                    u2 = nextu2;
                    fu2 = testcase1(u2);
                end
                final2 = u2;
                plot(0:600, final1, 'ro', 0:600, final2, 'b+', 0:600, uinitial, 'b-');
            else if (num == 4)
                    lambda = 0.8; 
                    tinit = 0; 
                    tend = 0.6;
                    gp = 41;
                    delx = 2/(gp-1);
                    x1 = -1;
                    for x=1:gp
                        r=x1+(x-1)*delx;
                        uinitial(x)=ux2(r);
                    end
                    delt = delx*lambda;
                    it = round((tend-tinit)/delt);
                    u = uinitial;
                    fu = testcase4(u);
                    n = length(uinitial);
                    for i=2:it
                        nextu(1) = 1/2*(u(n)+u(2)) + lambda/2*(fu(n) - fu(2));
                        for j=2:n-1
                            nextu(j) = 1/2*(u(j-1)+u(j+1)) + lambda/2*(fu(j-1) - fu(j+1));
                        end
                        nextu(n) = 1/2*(u(n-1)+u(1)) + lambda/2*(fu(n-1) - fu(1));
                        u = nextu;
                        fu = testcase4(u);
                    end
                    final = u;
                    plot(0:40, final, 'ro', 0:40, uinitial, 'b-');
                else if (num == 5)
                        lambda = 0.8; 
                        tinit = 0; 
                        tend = 0.3;
                        gp = 41;
                        delx = 2/(gp-1);
                        x1 = -1;
                        for x=1:gp
                            r=x1+(x-1)*delx;
                            uinitial(x)=ux5(r);
                        end
                        delt = delx*lambda;
                        it = round((tend-tinit)/delt);
                        u = uinitial;
                        fu = testcase4(u);
                        n = length(uinitial);
                        for i=2:it
                            nextu(1) = 1/2*(u(n)+u(2)) + lambda/2*(fu(n) - fu(2));
                            for j=2:n-1
                                nextu(j) = 1/2*(u(j-1)+u(j+1)) + lambda/2*(fu(j-1) - fu(j+1));
                            end
                            nextu(n) = 1/2*(u(n-1)+u(1)) + lambda/2*(fu(n-1) - fu(1));
                            u = nextu;
                            fu = testcase4(u);
                        end
                        final = u;
                        plot(0:40, final, 'ro', 0:40, uinitial, 'b-');
                    end
                end
            end
        end
    end
end