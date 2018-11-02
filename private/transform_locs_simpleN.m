function [transform,iAa,iBa]=transform_locs_simpleN(transform,channelref,locref,channeltarget,loctarget,p)
% locref , target: Mx2 or Mx3 array
sepscale=2;




if isfield(p,'Tfile') && exist(p.Tfile,'file')
    Tinitial=load(p.tfile);
    locT=Tinitial.transform2Reference(loctarget);
%     [loctT.x,loctT.y]=Tinitial.transformCoordinatesInv(loctarget.x(:),loctarget.y(:));
%     mirrorinfo=Tinitial.tinfo.mirror;
%     dx=0;dy=0;
else %all initial estimation:
    inforef=transform.info{channelref};
    infotarget=transform.info{channeltarget};
    locT(:,1)=loctarget(:,1)-infotarget.xrange(1);locT(:,2)=loctarget(:,2)-infotarget.yrange(1);
    locR(:,1)=locref(:,1)-inforef.xrange(1);locR(:,2)=locref(:,2)-inforef.yrange(1);
    
    %mirror if needed
    sref=[inforef.xrange(2)-inforef.xrange(1) inforef.yrange(2)-inforef.yrange(1)];
    for k=1:length(inforef.mirror)
        if inforef.mirror(k)>0
            locR(:,inforef.mirror(k))=sref(inforef.mirror(k))-locR(:,inforef.mirror(k));
        end
    end
    star=[infotarget.xrange(2)-infotarget.xrange(1) infotarget.yrange(2)-infotarget.yrange(1)];
    for k=1:length(inforef.mirror)
        if infotarget.mirror(k)>0
            locT(:,infotarget.mirror(k))=star(infotarget.mirror(k))-locT(:,infotarget.mirror(k));
        end
    end
    

    %determine approximate shift
    xr=1:1:sref(1);yr=1:sref(2);
    ht=histcounts2(locT(:,1),locT(:,2),xr,yr);
    hr=histcounts2(locR(:,1),locR(:,2),xr,yr);
    G=fftshift(ifft2(conj(fft2(hr)).*fft2(ht)));
    h=fspecial('gaussian',13,sepscale);
    Gf=filter2(h,G);
    [~ ,indmax]=max(Gf(:));
    [x0,y0]=ind2sub(size(Gf),indmax);
    dx0=x0-ceil(size(Gf,1)/2);
    dy0=y0-ceil(size(Gf,2)/2);
%     loctT.x=loctT.x-dx;
%     loctT.y=loctT.y-dy;
  

end
locRh.x=locR(:,1);locRh.y=locR(:,2);locRh.frame=ones(size(locRh.x));
locTh.x=locT(:,1);locTh.y=locT(:,2);locTh.frame=ones(size(locTh.x));
[iAa,iBa,na,nb,nseen]=matchlocsall(locRh,locTh,-dx0,-dy0,2*sepscale,1e5);

% transform=interfaces.LocTransform;
% t.type='polynomial';
% % t.type='affine';
% t.parameter=3;
transform.findTransform(channeltarget,locref(iAa,:),loctarget(iBa,:))



 tback=transform.transformToReference(channeltarget,loctarget(iBa,:));
dd=tback-locref(iAa,:);

%  figure(88);plot(tback(:,1),tback(:,2),'x',locref(:,1),locref(:,2),'o')  
%   figure(88);plot(locref.x,locref.y,'b.',loctT.x,loctT.y,'r+',loctT.x-dx0,loctT.y-dy0,'g.',loctargeti.x,loctargeti.y,'rx',xa,ya,'cx') 
   
if isfield(p,'ax')&& ~isempty(p.ax)
    axh=p.ax;
    plot(axh,dd(:,1),dd(:,2),'x')
    title(axh,[num2str(std(dd(:,1))) ', ' num2str(std(dd(:,2)))]);
end
% transform.tinfo.targetpos=targetpos;
% transform.tinfo.separator=separators;
% transform.tinfo.mirror=mirrorinfo;
% transform.tinfo.cam_pixelsize_nm=1;
% transform.tinfo.units='pixels';

end


function pos=reducepos(posin,df)
    z0=ceil(size(posin.x,1)/2);
    for l=size(posin.x,2):-1:1
        framerange=abs(posin.frame(:,l)-z0)<=df;
        x(:,l)=posin.x(framerange,l);y(:,l)=posin.y(framerange,l);z(:,l)=posin.z(framerange,l);frame(:,l)=posin.frame(framerange,l);
    end
    [~,indsort]=sort(frame(:));
    pos.x=x(indsort);pos.y=y(indsort);pos.z=z(indsort);pos.frame=frame(indsort);
    
end