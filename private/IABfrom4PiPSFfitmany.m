function out=IABfrom4PiPSFfitmany(PSF, phaseshift,frequency,roisizexy,numframes,PSFg)
% startpar=[dx dy dz phase norm] all vectors length 3 ;
% par= dxi
%      dyi
%      dzi
%      phii
%      Ni
%       %or rather normp for the four quadrants?
%      dxp
%      dyp
%      normp

% normf=par(10:12);
%shift PSF by dx,dy,dz

s=size(PSF);
rx=ceil((roisizexy+1)/2)+2;mp=ceil((s(1)+1)/2);mpz=ceil((s(3)+1)/2);rz=min(numframes,mpz-1);
PSFs=PSF(mp-rx:mp+rx,mp-rx:mp+rx,mpz-rz:mpz+rz,:,:);
Nistart=PSFg.globalnorm*ones(1,s(4));
dxstart=zeros(1,s(4));
dxpstart=-PSFg.dx;
dypstart=-PSFg.dy;
normpstart=PSFg.normf;
startpar=[dxstart(:);dxstart(:);dxstart(:);dxstart(:); Nistart(:); dxpstart(:);dypstart(:);normpstart(:)];
% startpar=[dxstart(:);dxstart(:);dxstart(:);dxstart(:); Nistart(:); dxpstart(:);dypstart(:);normpstart(:);phaseshift;frequency];
% startpar=[0 0 0 0 0 0 0 0 0 1 1 1];
% startpar=[0 0 0 0 0 0  1 1 1];
% fixpar=[phaseshift frequency];
% zshift0h=zshift0(2:4)-zshift0(1);
fixpar=[phaseshift frequency];

[PSFstart,mstart]=recoverPSF(startpar,PSFs,fixpar);

options=optimoptions('lsqcurvefit','Display','iter');
options.FinDiffRelStep=.01;
PSFc=PSFs(3:end-2,3:end-2,3:end-2,:,:);


fitpar=lsqcurvefit(@recoverPSF,startpar,PSFs,PSFc,[],[],options,fixpar);

[PSFrecovered,out]=recoverPSF(fitpar,PSF,fixpar);

out.dx=[0 fitpar(1:3)];
out.dy=[0 fitpar(4:6)];
out.dz=[0 zshift0h];
% out.dz=[0 fitpar(7:9)];
% out.normf=[1 fitpar(10:12)];
out.normf=[1 fitpar(7:9)];
out.frequency=frequency;
out.phaseshifts=[-pi phaseshift 0 phaseshift+pi];

end


function [PSFmSs,out]=recoverPSF(par,PSF,fixpar)
s=size(PSF);
numbeads=s(4);

dx=par(1:numbeads);
dy=par(numbeads+1:2*numbeads);
dz=par(2*numbeads+1:3*numbeads);
phi=par(3*numbeads+1:4*numbeads);
Ni=par(4*numbeads+1:5*numbeads);
off=5*numbeads;
dxp=par(off+1:off+4);
dyp=par(off+5:off+8);
normp=par(off+9:off+12);
% phaseshift=par(off+13);
% frequency=par(off+14);

dx(1)=0;
dy(1)=0;
dz(1)=0;
phi(1)=0;
dxp(1)=0;
dyp(1)=0;
normp(1)=1;

phaseshift=fixpar(1);
frequency=fixpar(2);

phaseshifts=[-pi phaseshift 0 phaseshift+pi];

%shift PSF and rescale
xn=1:size(PSF,1);yn=1:size(PSF,2);zn=1:size(PSF,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);
PSFS=zeros(size(PSF));
Ia=zeros(s(1),s(2),s(3),s(4));Aa=Ia;Ba=Ia;
for b=1:numbeads
    for c=1:4
        dxh=dx(b)+dxp(c);
        dyh=dy(b)+dyp(c);
        dzh=dz(b);
        PSFS(:,:,:,b,c)=interp3(PSF(:,:,:,b,c),Xq-dxh,Yq-dyh,Zq-dzh,'cubic',0);%/normf(k-1);
        [Ia(:,:,:,b),Aa(:,:,:,b),Ba(:,:,:,b)]=make4Pimodel(squeeze(PSFS(:,:,:,b,:)),phaseshifts+phi(b),frequency,normp*Ni(b));
    end
end

% calculate IAB
% [I,A,B]=make4Pimodel(PSFS,phaseshifts,frequency,[1 normf]);

%make PSF again
I=mean(Ia,4);A=mean(Aa,4);B=mean(Ba,4);
PSFm=makePSF(I,A,B,frequency, phaseshifts,  normp);

%shift back for fitting
PSFmS=zeros(size(PSF));
% PSFmS(:,:,:,1)=PSFm(:,:,:,1);
% for k=2:4
%     PSFmS(:,:,:,k)=interp3(PSFm(:,:,:,k),Xq+dx(k-1),Yq+dy(k-1),Zq+dz(k-1),'cubic',0);
% end

for b=1:numbeads
    for c=1:4
        dxh=dx(b)+dxp(c);
        dyh=dy(b)+dyp(c);
        dzh=dz(b);
        PSFmS(:,:,:,b,c)=interp3(PSFm(:,:,:,c),Xq+dxh,Yq+dyh,Zq+dzh,'cubic',0)*Ni(b);%/normf(k-1);
    end
end

PSFmSs=PSFmS(3:end-2,3:end-2,3:end-2,:,:);
% PSFm=PSFm(2:end-1,2:end-1,2:end-1,2:end-1);
out.PSF=PSFmS;
out.I=I;
out.A=A;
out.B=B;
end




function PSFo=makePSF(I,A,B,frequency, phaseshifts, normf)
s=size(I);
PSF=zeros(s(1)*s(2),s(3),4);
z=(1:s(3))-round(s(3)/2);
Ir=reshape(I,s(1)*s(2),s(3));
Br=reshape(B,s(1)*s(2),s(3));
Ar=reshape(A,s(1)*s(2),s(3));
for k=1:4
    PSF(:,:,k)=normf(k)*(Ir+Ar.*cos(2*frequency*z+phaseshifts(k))+Br.*sin(2*frequency*z+phaseshifts(k)));
end

PSFo=reshape(PSF,s(1),s(2),s(3),4);
end