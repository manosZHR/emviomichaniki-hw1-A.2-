function [dcdt]=odeA2(t,c,constants)

    ce=c(1); cs=c(2); ces=c(3); cp=c(4);
    k1=constants(1); k1r=constants(2); k2=constants(3);
    
    dcdt=zeros(4,1);
    
    dcdt(1)=(k1r+k2)*ces-k1*ce*cs;
    dcdt(2)=k1r*ces-k1*ce*cs;
    dcdt(3)=-dcdt(1);
    dcdt(4)=k2*ces;

end
