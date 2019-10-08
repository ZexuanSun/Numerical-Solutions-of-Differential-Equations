function lcrb(Dt,Tend,choice)
%左边中心右边节点
%10个点
%------------------------------------------
%step 1: set parameter: viacocity
%------------------------------------------
nu = 0.1;
%------------------------------------------
%step 2: define Spacial Domain and gride
%------------------------------------------
Xl = 0;
Xr = 1;
Ne = 9;%10个网格，
Nn = Ne + 1;
Dx = (Xr - Xl)/(Ne+1/2);
i = 1:Ne;
x = Dx/2 + Dx*(i-1);

time = 0;
new(2:Nn) = x.*(1-x);%最左边是鬼点,
                       %即new(1)
new(Nn+1) = 4*sin(6*time);
Crt = nu * Dt / Dx^2;
It = 1;
Nt = Tend/Dt;
while It <= Nt
    old = new;
    old(1) = old(2) - Dx * 10*sin(time);
 
    i = 2: Nn;
    new(i) = old(i) + Dt .* sin(2*pi.*x(i-1)) * sin(4*pi*time) + Crt*(old(i-1)- 2*old(i) + old(i+1));
    
    %new(1) = new(2) - Dx * 10*sin(time);
    time = time + Dt;
    new(Nn+1) = 4*sin(6*time);
    It = It + 1;
end
pnew = new(2:Nn+1);
x(Nn) = Xr;

%画10个点
[xe,ue,error] = matlab_sol(x,pnew,Dt,Tend);
if choice == 1
    %createfigure1(x,pnew,xe,ue)
     plot(x,pnew,'o--',xe,ue,'-') 
     legend('Num. Sol.','Pdepe. Sol.');
else
    
     plot(x,error,'o-')
     legend('Error');
end
sum(error)/size(x,2)
%plot(x,error,'o-')
% plot(x,pnew,'o--',xe,ue,'-') 
% legend('Num. Sol.','Pdepe. Sol.');
% %createfigure(x,pnew,xe,ue)
end