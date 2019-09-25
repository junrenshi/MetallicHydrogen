%% electron-test charge dielectric function with the local field correction

function y=q2epsUI_et(q,rs)

G=zeros(size(q));
y=zeros(size(q));
alp=(9*pi/4)^(1/3);
x=sqrt(rs);
b0=0.0621814;
b1=9.81379;
b2=2.82224;
b3=0.736411;
gam0=0.25+pi/(24*alp)*b0*(x^6*(2*x+4*b1*x^2+2*b1^2*x^3+4*b2*x^3+3*b1*b2*x^4+5*b3*x^4+4*b1*b3*x^5)/...
    (x^2+b1*x^3+b2*x^4+b3*x^5)^2+2*x^4*(1+b1*x)/(x^2+b1*x^3+b2*x^4+b3*x^5));
z=4*(rs/alp/pi)^(1/2);
g0=(z/besseli(1,z))^2/8;
A=0.029;
B=9.0/16*gam0-3.0/64*(1-g0)-16.0/15*A;
C=-3.0/4*gam0+9.0/16*(1-g0)-16.0/5*A;
ind=abs(q)<1e-6;
ind2=abs(q-2)<1e-8;
ind3=(abs(q-2)>=1e-8)&(abs(q)>=1e-6);
G(ind)=gam0*q(ind).^2;
% G(ind)=(2*B+8/3*A)*q(ind).^2;
G(ind2)=A.*q(ind2).^4+B*q(ind2).^2+C;
G(ind3)=A*q(ind3).^4+B*q(ind3).^2+C+(A*q(ind3).^4+(B+8/3*A)*q(ind3).^2-C).*...
    (4-q(ind3).^2)./(4*q(ind3)).*log(abs((2+q(ind3))./(2-q(ind3))));
y(ind)=q(ind).^2+4*(1-G(ind)).*rs./(alp*pi);
% y(ind)=q(ind).^2+4*rs./(alp*pi)./(1-4*rs./(alp*pi).*(2*B+8/3*A));
y(ind2)=q(ind2).^2+2*(1-G(ind2)).*rs./(alp*pi);
y(ind3)=q(ind3).^2+2*rs./(alp*pi).*(1-G(ind3)).*(1+(1-(q(ind3)/2).^2)./q(ind3).*...
    log(abs((1+q(ind3)/2)./(1-q(ind3)/2))));
