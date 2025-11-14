function [Vq,Vd,V0]=lipo_abc_to_qd0(Va,Vb,Vc,theta_e)
a=exp(1j*2*pi/3);
x_ab=(2/3)*(Va+a*Vb+a^2*Vc);
x_qd=x_ab.*exp(-1j*theta_e);
Vq=real(x_qd);
Vd=-imag(x_qd);
V0=(Va+Vb+Vc)/3;
end

