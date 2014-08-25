function [] = roe(num)
    if (num==1)    
        lambda = 0.8; 
        tinit = 0; 
        tend = 30;
        gp = 41;
        delx=2/(gp);
        x1 = -1;
        %uinitial = zeros(1:gp);
        for x=1:(gp)
            r=x1+(x-1)*delx;
            uinitial(x)=ux1(r);
        end
        delt = delx*lambda;
        it = round((tend-tinit)/delt);
        u = uinitial;
        fu = testf1(u);
        n = length(uinitial);
        for i=1:it
            nextu(1) = u(1) - lambda*(fu(1) - fu(n));
            for j=2:n-1
                nextu(j) = u(j) - lambda*(fu(j) - fu(j-1));
            end
            nextu(n) = u(n)- lambda*(fu(n) - fu(n-1));
            u = nextu;
            fu = testf1(u);
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
            fu = testf1(u);
            n = length(uinitial);
            for i=1:it
                nextu(1) = u(1) - lambda*(fu(1) - fu(n));
            for j=2:n-1
                nextu(j) = u(j) - lambda*(fu(j) - fu(j-1));
            end
            nextu(n) = u(n)- lambda*(fu(n) - fu(n-1));
                u = nextu;
                fu = testf1(u);
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
                fu1 = testf1(u1);
                n = length(uinitial);
                for i=2:iux1
                     nextu1(1) = u1(1) - lambda*(fu1(1) - fu1(n));
            for j=2:n-1
                nextu1(j) = u1(j) - lambda*(fu1(j) - fu1(j-1));
            end
            nextu1(n) = u1(n)- lambda*(fu1(n) - fu1(n-1));
                    u1 = nextu1;
                    fu1 = testf1(u1);
                end
                final1 = u1;
                delt = delx*lambda;
                u2 = uinitial;
                fu2 = testf1(u2);
                n = length(uinitial);
                for i=2:iux2
                    nextu2(1) = u2(1) - lambda*(fu2(1) - fu2(n));
            for j=2:n-1
                nextu2(j) = u2(j) - lambda*(fu2(j) - fu2(j-1));
            end
            nextu2(n) = u2(n)- lambda*(fu2(n) - fu2(n-1));
                    u2 = nextu2;
                    fu2 = testf1(u2);
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
                    fu = testf4(u);
                    n = length(uinitial);
                     for t=1:gp
                           r=x1+(t-1)*delx; 
                           if(t<=14) 
                            orig(t)=0.0;
                           else if( t<=26)
                            orig(t)= 1.67*(t-14)*delx;
                               else if(t<=33)
                                   orig(t) = 1.0;
                                       else
                                           orig(t)=0.0;
                                      end
                               end
                           end
                        end 
                    for i=2:it
                        nextu(1) = u(1) - lambda*min(0,(u(1)+u(2))/2)*(u(2)-u(1))-lambda*max(0,(u(1)+u(n))/2)*(u(1)-u(n));
                        for j=2:n-1
                            nextu(j) = u(j) - lambda*min(0,(u(j)+u(j+1))/2)*(u(j+1)-u(j))-lambda*max(0,(u(j)+u(j-1))/2)*(u(j)-u(j-1));
                        end
                        nextu(n) = u(n) - lambda*min(0,(u(n)+u(1))/2)*(u(1)-u(n))-lambda*max(0,(u(n)+u(n-1))/2)*(u(n)-u(n-1));
                        u = nextu;
                        fu = testf4(u);
                    end
                    final = u;
                    plot(0:40, final, 'ro', 0:40, orig, 'b-');
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
                        fu = testf4(u);
                        n = length(uinitial);
                        for i=2:it
                            nextu(1) = u(1) - lambda*min(0,(u(1)+u(2))/2)*(u(2)-u(1))-lambda*max(0,(u(1)+u(n))/2)*(u(1)-u(n));
                        for j=2:n-1
                            nextu(j) = u(j) - lambda*min(0,(u(j)+u(j+1))/2)*(u(j+1)-u(j))-lambda*max(0,(u(j)+u(j-1))/2)*(u(j)-u(j-1));
                        end
                        nextu(n) = u(n) - lambda*min(0,(u(n)+u(1))/2)*(u(1)-u(n))-lambda*max(0,(u(n)+u(n-1))/2)*(u(n)-u(n-1));
                            u = nextu;
                            fu = testf4(u);
                        end
                        final = u;
                        for t=1:gp
                           r=x1+(t-1)*delx; 
                           if(t<=8) 
                            orig(t)=-1.0;
                           else if( t<=20)
                            orig(t)= -1.0 + 2/12*(t-8);
                               else if(t<=27)
                                   orig(t) = 1;
                                       else
                                           orig(t)=-1.0;
                                      end
                               end
                           end
                        end
                        plot(0:40, final, 'ro', 0:40,orig, 'b-');
                    end
                end
            end
        end
    end
end
