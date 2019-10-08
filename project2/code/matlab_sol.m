function [x,u,error] = matlab_sol(y,py,Dt,Tend)
m = 0;
 num = size(y,2);
if size(y,2) == 10
    y(11) = 0;
    y = sort(y);
end
x = [];
for i = 1:10
    x = [x linspace(y(i), y(i+1),11)];
end
x = [x linspace(y(11),1,11)];
x = unique(x);
%x = y;
%x = linspace(0,1,101);
% x = sort([x y]);
% x = unique(x);
%x = sort([x,y]);
spac = Tend/Dt;
t = linspace(0,Tend,spac);

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);
u = u(end,:);

if num == 10
    for i = 1:10
        %error(i) = u(i*10+1) - py(i);
        error(i) = abs(u(i*10+1) - py(i));
    end
else
    for i = 1:11
        %error(i) = u((i-1)*10+1) - py(i);
        error(i) = abs(u((i-1)*10+1) - py(i));
    end
end


% % A surface plot is often a good way to study a solution.
% surf(x,t,u) 
% title('Numerical solution computed with 20 mesh points.')
% xlabel('Distance x')
% ylabel('Time t')

% A solution profile can also be illuminating.
% figure
% plot(x,u(end,:))
% title('Solution at t = 0.9')
% xlabel('Distance x')
% ylabel('u(x,2)')
end
% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
c = 10;
f = DuDx;
s = 10 * sin(2*pi*x)*sin(4*pi*t);
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = x*(1-x);
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = -10*sin(t);
ql = 1;
pr = ur - 4*sin(6*t);
qr = 0;
end