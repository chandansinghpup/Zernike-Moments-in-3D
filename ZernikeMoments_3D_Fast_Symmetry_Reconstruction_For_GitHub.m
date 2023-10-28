
% Computation of 3D Zernike moments of a function f(x,y,z) upto order pmax
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  Authors: Dr. Chandan Singh, Professor(Retd.),
 %           Dr. Sukhjeet Kaur Ranade, Professor,
 %           Rama Rani,
 %           Rupneet Kaur,
 %           Department of Computer Science,   
 %          Punjabi University, Patiala, Punjab, 147002, India.
 %          email: chandan.sp@gmail.com
 %  Date: 20-04-2023
 %
 %  References : An-Wen Deng, Chin-Yen Gwo, A Stable Algorithm Computing
 %               High Order 3D Zernike Moments and Shape Reconstruction, 
 %               Proceedings of the 2020 4th International Conference on 
 %               Digital Signal Processing 2020, pp. 38-42.
 %  
 %               An-Wen Deng, Chin-Yen Gwo, Efficient and Stable Voxel-Based 
 %               Algorithm for Computing 3D Zernike Moments and Shape 
 %               Reconstruction2021.
 %  
 %               Jérôme Houdayer, Patrice Koehl,Stable Evaluation of 3D Zernike
 %               Moments for Surface Meshes, Algorithms, 15(11), p.406.
 %               
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;  clc;              
tstart=tic;
pmax=50;  option_circle=1;  % 1-outer, 2-inner
% Compute factorials upto pmax+1 and save.
pmaxp1=pmax+1;  step=2;  % order of 3D ZMs.

% A synthetic image.
M=128; 
% M=8
N=M;  T=M;  denom=(M-1)*(N-1)*(T-1);
f=zeros(M,N,T);
for x=1:M
    for y=1:N
        for z=1:T
            f(x,y,z)=255.0*(x-1)*(y-1)*(z-1)/denom;
        end
    end
end
[M,N,T]=size(f);  
fprintf(1, 'Image size: M=%d, N=%d, T=%d\n', M, N,T);
f=double(f);
T2=floor(T/2);  Mp1=M+1;  Np1=N+1;  Tp1=T+1;

if(mod(M,2)==1)
    error('myApp:argChk', ['Size of image, M=%d is not even. \n8-way '...
        'symmetry works for even-sized image'], M);
end
if(~(M==N && N==T))
    error('myApp:argChk', ['Size of image [%d, %d, %d], is not the same in '...
        ' all dimensions. \n8-way '...
        'symmetry works for even-sized image with equal size in \n'...
        'all dimensions'], M, N, T);
end
avg=mean(f,'all');  fmin=min(f,[],'all');  fmax=max(f,[],'all');
fprintf(1, 'M=%d, N=%d, T=%d, pmax=%d, average=%f, fmin=%f, fmax=%f\n', ...
    M, N, T, pmax, avg, fmin, fmax);  

% Compute the coefficients of Legendre polynomials(LP) of
% degree l, and repetition m, recursively.
LP=zeros(pmaxp1, pmaxp1);
L=pmax*(pmax+3)/2;
C1=zeros(L,1);  C2=zeros(L,1);  C3=zeros(L,1);  C4=zeros(L,1);
C5=zeros(L,1);
LP(1,1)=1.0;  % for l=0, LP(0,0)=1.0
indx=0;
for l=1:pmax  % for m==l
    indx=indx+1;
    C3(indx)=sqrt((2*l+1)/(2*l));
end
for l=1:pmax  % for m>1 and 0<=m<=l-1, i.e., m=0,1,2,...,l-2
    for m=0:l-1
        indx=indx+1;
        if(indx>L)
            error('myApp:argChk', ['indx=%d, which is more than maximum '
                'allowed value, L=%d'], indx, L);
        end
        temp1=sqrt((2*l+1)/((l+m)*(l-m)));
        C1(indx)=temp1*sqrt(2*l-1);
        temp2=(l+m-1)*(l-m-1)/(2*l-3);
        if(temp2<0.0)
            error('myApp:argChk', ['squar root of temp2 '...
                'which is -ve, temp2=%f\n'], temp2);
        end
        C2(indx)=temp1*sqrt(temp2);
    end  % next m
end  % next l

% Compute coefficients of Zernike radial functions and save.
if(mod(pmax,2)==1)
    L=(pmax+1)*(pmax+3) /4-2;  
else L=(pmax+2)*(pmax+2)/4-2;  
end

K1=zeros(L,1);  K2=zeros(L,1);  K3=zeros(L,1);  
indx=0;  
for q=0:pmax
    qp1=q+1;
    for p=4+q:step:pmax
        indx=indx+1; 
        if(indx>L)
            error('myApp:argChk', ['indx=%d, which is more than maximum '
                'allowed value, L=%d'], indx, L);
        end
        pp1=p+1;
        k0=double((p-q)*(p+q+1)*(2*p-3));
        k1=(2*p-1)*(2*p+1)*(2*p-3);
        k2=-0.50*((2*p-1)*(4*q*(q+1)+1)+k1);
        k3=-(p-q-2)*(p+q-1)*(2*p+1);
        K1(indx)=k1/k0;  K2(indx)=k2/k0; K3(indx)=k3/k0;
    end  % next p
end  % next q

if(option_circle==1) D=sqrt(double(M*M+N*N+T*T));
elseif(option_circle==2) D=double(min(M,N,T));
end
RAD=floor(M/2);  RADp1=RAD+1;
coeff=2.0/(pi*D*D*D);  % if M=N=T, then D*D*D=3*sqrt(3)*M*M*M.
% No. of moments
npqm=0;
for p=0:pmax
    for q=p:-step:0
        for m=0:q
            npqm=npqm+1;
        end
    end
end
NPQM=npqm;
fprintf(1,'\nNo. of moments, NPQM=%d\n', NPQM);
ZR=zeros(NPQM,1);  ZI=zeros(NPQM,1);  ZM=zeros(NPQM,1);  
rpower=zeros(pmaxp1,1);  cosphiqm=zeros(pmaxp1,1);  sinphiqm=zeros(pmaxp1,1);  
for X=1:RAD
    if(mod(X-1,10)==0) fprintf(1,'Computing ZMs for X=%d to X=%d\n', X, X+10);  end
    xf=(2.0*X-1)/D;  xf2=xf*xf;
    X1=X+RAD;  X2=RADp1-X; 
    for Y=1:X  % <=Y<=X
        Y1=Y+RAD;  Y2=RADp1-Y; 
        yf=(2.0*Y-1)/D;  yf2=yf*yf;
        rxy=sqrt(xf2+yf2);
        if(rxy>1.0) continue;  end
        if(rxy<eps)  cosphi=1.0;  sinphi=0.0;  % to avoid division by zero.
        else cosphi=xf/rxy;  sinphi=yf/rxy;
        end
        
        % Recursive computation of sine(phi) and cosine(phi) functions.
        cosphiqm(1)=1.0;  sinphiqm(1)=0.0;   % for q=0.0.
        for m=1:pmax
            mp1=m+1;
            cosphiqm(mp1)=cosphiqm(m)*cosphi-sinphiqm(m)*sinphi;
            sinphiqm(mp1)=sinphiqm(m)*cosphi+cosphiqm(m)*sinphi;
        end  % next m
        for Z=1:RAD
            Z1=Z+RAD;  Z2=RADp1-Z;
            zf=(2.0*Z-1)/D;  zf2=zf*zf;
            r=sqrt(xf2+yf2+zf2);
            if(r>1.0) continue;  end
            
            if(r<eps)  cost=1.0;  sint=0.0;  % to avoid division by zero.
            else cost=zf/r;  sint=rxy/r;
            end
            
            f1=f(X1,Y1,Z1);  f2=f(X2,Y1,Z1);  f3=f(X2,Y2,Z1);  f4=f(X1,Y2,Z1);
            f5=f(Y1,X1,Z1);  f6=f(Y1,X2,Z1); f7=f(Y2,X2,Z1); f8=f(Y2,X1,Z1);
            f9=f(X1,Y1,Z2);  f10=f(X2,Y1,Z2);  f11=f(X2,Y2,Z2);  f12=f(X1,Y2,Z2);
            f13=f(Y1,X1,Z2); f14=f(Y1,X2,Z2); f15=f(Y2,X2,Z2); f16=f(Y2,X1,Z2);
            
            F1=f1+f2+f3+f4;      F2=f9+f10+f11+f12;    F3=f5+f6+f7+f8;
            F4=f13+f14+f15+f16;
            F5=-f1+f2-f3+f4;     F6=-f9+f10-f11+f12;   F7=+f5-f6+f7-f8;
            F8=f13-f14+f15-f16;
            F9=f1-f2-f3+f4;      F10=-f9+f10+f11-f12;  F11=f5+f6-f7-f8;  
            F12=-f13-f14+f15+f16;
            F13=-f1-f2+f3+f4;    F14=f9+f10-f11-f12;   F15=-f5+f6+f7-f8;  
            F16=f13-f14-f15+f16;
            
            % Recursive computue r**p
            mult=r;
            rpower(1)=1.0;
            for p=1:pmax
                rpower(p+1)=mult;
                mult=mult*r;
            end  % next p
            
            % Computaion of moments
            
            % Compute Legendre polynomials(LP) of degree l recursively.
            indx=0;
            for l=1:pmax  % for m==l
                indx=indx+1;
                LP(l+1,l+1)=-C3(indx)*sint*LP(l-1+1,l-1+1);% LP(l,l)=-C3*sin(theta)
                %*LP(l-1,l-1)
            end
            for l=1:pmax  % for m>1 and 0<=m<=l-1, i.e., m=0,1,2,...,l-2
                for m=0:l-1
                    indx=indx+1;
                    if(l==1)
                        LP(l+1,m+1)=C1(indx)*cost*LP(l-1+1,m+1);  % LP(-1,0)=0.0
                    else
                        LP(l+1,m+1)=C1(indx)*cost*LP(l-1+1,m+1)-C2(indx)*LP(l-2+1,m+1);
                    end
                end  % next m
            end  % next l
            
            % Derive Zernike radial polynomials ZRP(p,q) of order p and
            % repeation q.
            ZRP=zeros(L,1);
            %fxyz=f(x,y,z);
            for p=0:pmax
                pp1=p+1;
                ZRP(pp1,pp1)=rpower(pp1);
                if(p<=1) continue;
                else ZRP(pp1,pp1-2)=...
                        0.5*((2*p+1)*rpower(pp1)-(2*p-1)*rpower(pp1-2));
                end
            end  % next p
            indx=0;
            for q=0:pmax
                qp1=q+1;
                for p=4+q:step:pmax
                    indx=indx+1;
                    pp1=p+1;
                    ZRP(pp1,qp1)=(K1(indx)*r*r+K2(indx))*ZRP(pp1-2,qp1)+...
                        K3(indx)*ZRP(pp1-4,qp1);
                end  % next p
            end  % next q
            
            % Compute 3D ZMs.
            npqm=0;  sign=-1;
            for p=0:pmax
                pp1=p+1;  
                sign=-sign;  % Since p starts with 0, first time it is +ve
                %if(mod(p,2)==0) sign=1;  else sign=-1;  end
                for q=p:-step:0
                    %sign=(-1)^q;
                    qp1=q+1;
                    value=ZRP(pp1,qp1);  
                    for m=0:q
                        npqm=npqm+1;
                        mp1=m+1;  rem=mod(m,4);
                        value1=value*LP(qp1,mp1);
                        if(X==Y)  % 4-way symmetry only
                            switch rem
                                case 0
                                    FR=F1+sign*F2;
                                    FI=F5+sign*F6;
                                    smR=FR*cosphiqm(mp1);
                                    smI=FI*sinphiqm(mp1);
                                case 1
                                    FR=F9+sign*F10;
                                    FI=F13+sign*F14;
                                    smR=FR*cosphiqm(mp1);
                                    smI=FI*sinphiqm(mp1);
                                case 2  % current 4-way modified from 8-way
                                    FR=F1+sign*F2;
                                    FI=F5+sign*F6;
                                    smR=FR*cosphiqm(mp1);
                                    smI=FI*sinphiqm(mp1);
                                case 3
                                    FR=F9+sign*F10;
                                    FI=F13+sign*F14;
                                    smR=FR*cosphiqm(mp1);
                                    smI=FI*sinphiqm(mp1);
                                otherwise
                                    error('myApp:argChk', ['Value of m=%d, '...
                                        'which is wrong'], m);
                            end  % of switch
                        else  % symmetry for 8 points
                            switch rem
                                case 0
                                    FR=F1+F3+sign*(F2+F4);
                                    FI=F5+F7+sign*(F6+F8);
                                    smR=FR*cosphiqm(mp1);
                                    smI=FI*sinphiqm(mp1);
                                case 1
                                    FR1=F9+sign*F10;
                                    FR2=F11+sign*F12;
                                    FI1=F13+sign*F14;
                                    FI2=F15+sign*F16;
                                    smR=(FR1*cosphiqm(mp1)+FR2*sinphiqm(mp1));
                                    smI=(FI1*sinphiqm(mp1)+FI2*cosphiqm(mp1));
                                case 2
                                    FR=F1-F3+sign*(F2-F4);
                                    FI=F5-F7+sign*(F6-F8);
                                    smR=FR*cosphiqm(mp1);
                                    smI=FI*sinphiqm(mp1);
                                case 3
                                    FR1=F9+sign*F10;
                                    FR2=-(F11+sign*F12);
                                    FI1=F13+sign*F14;
                                    FI2=-(F15+sign*F16);
                                    smR=(FR1*cosphiqm(mp1)+FR2*sinphiqm(mp1));
                                    smI=(FI1*sinphiqm(mp1)+FI2*cosphiqm(mp1));
                                otherwise
                                    error('myApp:argChk', ['Value of m=%d, '...
                                        'which is wrong'], m);
                            end  % of switch
                        end  % if X=Y
                        ZR(npqm)=ZR(npqm)+value1*smR;
                        ZI(npqm)=ZI(npqm)+value1*smI;
                    end  % next m
                end  % next q
            end  % next p
        end  % next z
    end  % next y
end  % next X

npqm=0;
for p=0:pmax
    pp1=p+1;  coeff1=coeff*(2*p+3);
    for q=p:-step:0
        qp1=q+1;
        for m=0:q
            mp1=m+1;  
            npqm=npqm+1;
            if(npqm>NPQM)
                error('myApp:argChk', ['No. of moments, '...
                    'npqm=%d is more than the actual value, '...
                    'NPQM=%d'], npqm, NPQM);
            end
            ZRpq=coeff1*ZR(npqm);  ZR(npqm)=ZRpq;
            ZIpq=coeff1*ZI(npqm);  ZI(npqm)=ZIpq;
            ZMpq=sqrt(ZRpq*ZRpq+ZIpq*ZIpq);  ZM(npqm)=ZMpq;
            if p==10
                fprintf(1,'displaying ZMs for pmax=%d\n', 10);
                fprintf(1, 'p=%d, q=%d, m=%d, ZRpq=%f, ZIpq=%f, ZMpq=%f\n',...
                p, q, m, ZRpq, ZIpq, ZMpq);
            end
        end  % next m
    end  % next q
end  % next p
fprintf(1, '\nTotal time taken for ZMs Computation=%f(sec)\n ', toc(tstart));

% 3D Image reconstruction.
tstart1=tic;  fR=zeros(M,N,T); foriginal=zeros(M,N); fRxy=zeros(M,N); 

for X=1:RAD
    if(mod(X-1,10)==0) fprintf(1,'Recostructing image for X=%d to X=%d\n', X, X+10);  end
    xf=(2.0*X-1)/D;  xf2=xf*xf;
    X1=X+RAD;  X2=RADp1-X; 
    for Y=1:X  % <=Y<=X
        Y1=Y+RAD;  Y2=RADp1-Y; 
        yf=(2.0*Y-1)/D;  yf2=yf*yf;
        rxy=sqrt(xf2+yf2);
        if(rxy>1.0) continue;  end
        if(rxy<eps)  cosphi=1.0;  sinphi=0.0;  % to avoid division by zero.
        else cosphi=xf/rxy;  sinphi=yf/rxy;
        end
        
        % Recursive computation of sine(phi) and cosine(phi) functions.
        cosphiqm(1)=1.0;  sinphiqm(1)=0.0;   % for q=0.0.
        for m=1:pmax
            mp1=m+1;
            cosphiqm(mp1)=cosphiqm(m)*cosphi-sinphiqm(m)*sinphi;
            sinphiqm(mp1)=sinphiqm(m)*cosphi+cosphiqm(m)*sinphi;
        end  % next m
        for Z=1:RAD
            Z1=Z+RAD;  Z2=RADp1-Z;
            zf=(2.0*Z-1)/D;  zf2=zf*zf;
            r=sqrt(xf2+yf2+zf2);
            if(r>1.0) continue;  end
            if(r<eps)  cost=1.0;  sint=0.0;  % to avoid division by zero.
            else cost=zf/r;  sint=rxy/r;
            end
            % Recursive computue r**p
            mult=r;
            rpower(1)=1.0;
            for p=1:pmax
                rpower(p+1)=mult;
                mult=mult*r;
            end  % next p
            
            indx=0;
            for l=1:pmax  % for m==l
                indx=indx+1;
                LP(l+1,l+1)=-C3(indx)*sint*LP(l-1+1,l-1+1);% LP(l,l)=-C3*sin(theta)
                %*LP(l-1,l-1)
            end
            for l=1:pmax  % for m>1 and 0<=m<=l-1, i.e., m=0,1,2,...,l-2
                for m=0:l-1
                    indx=indx+1;
                    if(l==1)
                        LP(l+1,m+1)=C1(indx)*cost*LP(l-1+1,m+1);%LP(-1,0)=0.0
                    else
                        LP(l+1,m+1)=C1(indx)*cost*LP(l-1+1,m+1)-...
                            C2(indx)*LP(l-2+1,m+1);
                    end
                end  % next m
            end  % next l
            
            % Derive Zernike radial polynomials ZRP(p,q) of order p and
            % repeation q.
            ZRP=zeros(pmax+1,pmax+1);
            for p=0:pmax
                pp1=p+1;
                ZRP(pp1,pp1)=rpower(pp1);
                if(p<=1) continue;
                else ZRP(pp1,pp1-2)=...
                        0.5*((2*p+1)*rpower(pp1)-(2*p-1)*rpower(pp1-2));
                end
            end  % next p
            indx=0;
            for q=0:pmax
                qp1=q+1;
                for p=4+q:step:pmax
                    pp1=p+1;  indx=indx+1;
                    ZRP(pp1,qp1)=(K1(indx)*r*r+K2(indx))*ZRP(pp1-2,qp1)+...
                        K3(indx)*ZRP(pp1-4,qp1);
                end  % next q
            end  % next p
            
            % Image reconstruction.
            f1=0.0;  f2=0.0;  f3=0.0;  f4=0.0;  f5=0.0;  f6=0.0;  f7=0.0;  
            f8=0.0;  f9=0.0; f10=0.0;  f11=0.0; f12=0.0;  f13=0.0;  
            f14=0.0;  f15=0.0;    f16=0.0; 
            npqm=0;  signq=-1;
            for p=0:pmax
                pp1=p+1;
                signq=-signq;  % Since p starts with 0, first time it is +ve
                for q=p:-step:0
                    qp1=q+1;  
                    value=ZRP(pp1,qp1);
                    sum_m=0.0;  signm=-1;   
                    for m=0:q
                        npqm=npqm+1;
                        signm=-signm;  % Since m starts with 0, first time it is +ve
                        signqm=signq*signm;
                        if(npqm>NPQM)
                            error('myApp:argChk', ['No. of moments, '...
                                'npqm=%d is more than the actual value, '...
                                'NPQM=%d'], npqm, NPQM);
                        end
                        mp1=m+1;  
                        cosphim=cosphiqm(mp1);  sinphim=sinphiqm(mp1); 
                        ZRpq=ZR(npqm);
                        ZIpq=ZI(npqm);
                        value1=value*LP(qp1,mp1);
                        F1=ZRpq*cosphim;  F2=ZIpq*sinphim;  
                        F3=ZRpq*sinphim;  F4=ZIpq*cosphim;  
                        S1=value1*(F1+F2);  S2=value1*(F1-F2);  
                        S3=value1*(F3+F4);  S4=value1*(F3-F4); 
                        % Both for 4-way and 8-way symmetry
                        if(m==0)
                            f1=f1+S2;  f2=f2+S1;
                            f3=f3+S2;  f4=f4+S1;
                            f9=f9+signq*S2;  f10=f10+signq*S1;
                            f11=f11+signq*S2;  f12=f12+signq*S1;
                        else
                            f1=f1+2.0*S2;  f2=f2+2.0*signm*S1;
                            f3=f3+2.0*signm*S2;  f4=f4+2.0*S1;
                            f9=f9+2.0*signqm*S2;  f10=f10+2.0*signq*S1;
                            f11=f11+2.0*signq*S2;  f12=f12+2.0*signqm*S1;
                        end
                        if(~(X==Y))  % for 8-way symmetry
                            switch mod(m,4)
                                case 0
                                    if(m==0)
                                        f5=f5+S1;
                                        f6=f6+S2;
                                        f7=f7+S1;  
                                        f8=f8+S2;
                                        f13=f13+signq*S1;
                                        f14=f14+signq*S2;
                                        f15=f15+signq*S1;
                                        f16=f16+signq*S2;
                                    else
                                        f5=f5+2.0*S1;
                                        f6=f6+2.0*S2;
                                        f7=f7+2.0*S1;  
                                        f8=f8+2.0*S2;
                                        f13=f13+signq*2.0*S1;
                                        f14=f14+signq*2.0*S2;
                                        f15=f15+signq*2.0*S1;
                                        f16=f16+signq*2.0*S2;
                                    end
                                case 1
                                    f5=f5+2.0*S4;  
                                    f6=f6+2.0*S3;
                                    f7=f7-2.0*S4;  
                                    f8=f8-2.0*S3;
                                    f13=f13-signq*2.0*S4;
                                    f14=f14-signq*2.0*S3;
                                    f15=f15+signq*2.0*S4;
                                    f16=f16+signq*2.0*S3;
                                case 2
                                    f5=f5-2.0*S1;  
                                    f6=f6-2.0*S2;
                                    f7=f7-2.0*S1;  
                                    f8=f8-2.0*S2;
                                    f13=f13-signq*2.0*S1;
                                    f14=f14-signq*2.0*S2;
                                    f15=f15-signq*2.0*S1;
                                    f16=f16-signq*2.0*S2;
                                case 3
                                    f5=f5-2.0*S4;  
                                    f6=f6-2.0*S3;
                                    f7=f7+2.0*S4;  
                                    f8=f8+2.0*S3;
                                    f13=f13+signq*2.0*S4;
                                    f14=f14+signq*2.0*S3;
                                    f15=f15-signq*2.0*S4;
                                    f16=f16-signq*2.0*S3;
                                otherwise
                                    error('myApp:argChk', ['Value of m=%d, '...
                                        'which is wrong'], m);
                            end  % of switch
                        end  % of if(~(X==Y))
                    end  % next m
                end  % next q
            end  % next p
            fR(X1,Y1,Z1)=f1;  fR(X2,Y1,Z1)=f2;  fR(X2,Y2,Z1)=f3;  fR(X1,Y2,Z1)=f4;
            fR(X1,Y1,Z2)=f9;  fR(X2,Y1,Z2)=f10;  fR(X2,Y2,Z2)=f11;  fR(X1,Y2,Z2)=f12;
            if(~(X==Y))
                fR(Y1,X1,Z1)=f5; fR(Y1,X2,Z1)=f6; fR(Y2,X2,Z1)=f7; 
                fR(Y2,X1,Z1)=f8;  fR(Y1,X1,Z2)=f13; fR(Y1,X2,Z2)=f14; 
                fR(Y2,X2,Z2)=f15; fR(Y2,X1,Z2)=f16;
            end  % if(~(X==Y))
        end  % next Z
    end  % next Y
end  % next X
foriginalT2(:,:)=f(:,:,T2);
foriginalT2p1(:,:)=f(:,:,T2+1);
fRxyT2(:,:)=fR(:,:,T2);
fRxyT2p1(:,:)=fR(:,:,T2+1);
NMSE=sum((f-fR).^2, 'all')/sum((f).^2, 'all');
fprintf(1,'\nNormalized Mean Square Error(NMSE)of 3D image=%f\n', NMSE);
fprintf(1, '\nTotal time taken for ZMs reconstruction=%f(sec)\n\n',...
    toc(tstart1));
fprintf(1, ['\nTotal time taken for ZMs Computation and '...
    'reconstruction=%f(sec)\n'], toc(tstart));
figure;
subplot(2,2,1), imshow(foriginalT2, []), title('Original middle 2D slice');
subplot(2,2,2), imshow(foriginalT2p1, []), title('Original middle+1 2D slice');
subplot(2,2,3), imshow(fRxyT2, []), title('Reconstructed middle 2D slice');
subplot(2,2,4), imshow(fRxyT2p1, []), title('Reconstructed middle+1 2D slice');


 
% Generate the isosurface for the input image (iso1)
iso1 = isosurface(f);
% Generate the isosurface for the reconstructed image (iso_reconstruct)
% Convert negative pixel values to zero and normalize to the range [0, 255]
min_fR = min(fR, [], 'all');
max_fR = max(fR, [], 'all');
range = max_fR - min_fR;
fRR = 255.0 * (fR - min_fR) / range;
% Generate the isosurface for the reconstructed image
iso_reconstruct = isosurface(fRR);
% Create a figure to display both plots side by side
figure;
% First subplot - 3D Volume Visualization of Brain Image (Input image)
subplot(1, 2, 1); % Create a subplot with 1 row and 2 columns, select the first plot (input image)
p1 = patch(iso1);
set(p1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); % Customize colors
daspect([1 1 1]); % Equal aspect ratio
view(3); % 3D view
axis tight;
camlight; lighting gouraud;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Input image');
grid on;
% Second subplot - 3D Zernike Volume Visualization of Brain Image (Reconstructed image)
subplot(1, 2, 2); % Select the second plot (reconstructed image)
p2 = patch(iso_reconstruct);
set(p2, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); % Customize colors
daspect([1 1 1]); % Equal aspect ratio
view(3); % 3D view
axis tight;
camlight; lighting gouraud;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Reconstructed image');
grid on;

