close all;

set(0,'defaultlinelinewidth',1.5)
%Constants
h = 6.626e-34;
c = 2.998e8;
h_cut = 1.055e-34;
m0 = 9.109e-31;
e_const = 1.602e-19;

%User inputs
U_eV=1;
a=3e-10;
b=4e-10;

%Derived values
U=U_eV*e_const;
ap0 = sqrt(2*m0*U/(h_cut^2));

f=@(g) (1-2*g)./(2*sqrt(g.*(g-1))).*sin(a*ap0*sqrt(g))...
    .*sin(b*ap0*sqrt(g-1))...
    +cos(a*ap0*sqrt(g)).*cos(b*ap0*sqrt(g-1));


g=linspace(.1, 10, 1e5);
fg=f(g);
g(isnan(g))=1;
plot(g,fg)
hold on
plot([g(1) g(end)], [1, 1], 'r--')
plot([g(1) g(end)], [-1, -1], 'g--')
ylim([min(fg)-.5, 3])
xlabel('\zeta (=E/U_o) \rightarrow');
ylabel('f(\zeta) (=RHS) \rightarrow');
title(['Plot of the RHS of the equation f(\zeta) vs. \zeta; '...
    'for U=' num2str(U_eV), ' eV; a=' num2str(a*1e10) char(197)...
    ' and b='  num2str(b*1e10) char(197)])
grid on

flg=abs(fg)<=1;
figure
h1=gca;
hold on
xlabel('Crystal momentum, k(radian/meter) \rightarrow');
ylabel('Energy, E (eV) \rightarrow');
title(['Reduced zone representation of the E-k relationship '...
    'for U=' num2str(U_eV), ' eV; a=' num2str(a*1e10) char(197)...
    ' and b='  num2str(b*1e10) char(197)])
xticks([-pi -pi/2 0 pi/2 pi]/(a+b));
xticklabels({'-\pi/(a+b)','-\pi/(2(a+b))','0','\pi/(2(a+b))','\pi/(a+b)'})
grid on

figure
h2=gca;
hold on
xlabel('Crystal momentum, k(radian/meter) \rightarrow');
ylabel('Energy, E (eV) \rightarrow');
title(['Extended zone representation of the E-k relationship '...
    'for U=' num2str(U_eV), ' eV; a=' num2str(a*1e10) char(197)...
    ' and b='  num2str(b*1e10) char(197)])
xticks([-6*pi -5*pi -4*pi -3*pi -2*pi -pi 0 pi 2*pi 3*pi 4*pi 5*pi 6*pi]/(a+b))
xticklabels({'-6\pi/(a+b)','-5\pi/(a+b)','-4\pi/(a+b)' ...
    '-3\pi/(a+b)','-2\pi/(a+b)','-\pi/(a+b)','0','\pi/(a+b)',... 
    '2\pi/(a+b)','3\pi/(a+b)'...
    '4\pi/(a+b)','5\pi/(a+b)','6\pi/(a+b)'})
xtickangle(45)
grid on

prd=pi/(a+b);
plst=1;
k=1;
while ~isempty(flg) && k<6
    pos=find(flg);
    if isempty(pos)
        break
    end
    pfst=plst+pos(1)-1;
    flg=flg(pos(1):end);
    pos=find(~flg);
    if isempty(pos)
        break
    end
    plst=pfst+pos(1)-1;
    flg=flg(pos(1):end);
    
    kv=acos(fg(pfst:plst-1))/(a+b);
    ev=g(pfst:plst-1)*U_eV;
    if mod(k,2)
        plot(h1,[-fliplr(kv), kv], [fliplr(ev),ev], 'b');
        if k==1
            plot(h2,[-fliplr(kv), kv], [fliplr(ev),ev], 'b');
        else
            plot(h2,kv+prd*(k-1), ev, 'b');
            plot(h2,-fliplr(kv)-prd*(k-1), fliplr(ev), 'b');
        end
    else
        plot(h1, [kv, -fliplr(kv)], [ev,fliplr(ev)], 'b');
        plot(h2,kv-prd*k, ev, 'b');
        plot(h2,-fliplr(kv)+prd*k, fliplr(ev), 'b')
    end
    k=k+1;
end
