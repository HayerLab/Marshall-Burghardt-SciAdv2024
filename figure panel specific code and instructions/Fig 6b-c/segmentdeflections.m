function mask=segmentdeflections(mask,nucr,cellsizeratio,debrisarea)
%%% only consider big nuclei %%%%%
bigthresh=round(cellsizeratio*pi*nucr^2);  %1.33*pi*nucr^2 is the avg cell size
onlybig=bwareaopen(mask,bigthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucedges=bwmorph(onlybig,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(mask));
%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    coorset=[c,r];                      %adjust to x-y convention
    [order,status]=orderperimeter(coorset);
    if status==0                    %error, skip segmentation for this cell
        continue
    end
    orderedset=coorset(order,:);
    bordermask=splitdeflections(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
mask=bwareaopen(mask,debrisarea); %bin1:0.75 bin2:0.5
end

function [bordermask,bridgeflag]=splitdeflections(orderedset,bordermask,nucr)
bridgeimage=0;
bridgeflag=0; %returned as 1 if deflections are bridged
nucr=round(nucr/4)*4;
%minperi=50;  %minimum length of perimeter to consider breaking

%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
offsetshort=nucr/4;
if offsetshort==1
    offsetshort=2;
end
offsetlong=2*offsetshort;
gradientoffsetshort=offsetshort;
gradientoffsetlong=offsetlong;
shortgradthresh = pi/6;
longgradthresh = pi/6;

%%% calculate deflections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numpix=size(orderedset,1);
%%%%%%% short steps %%%%%%%%%%%%%%%%%%
orderedsetoffsetshort=[orderedset(offsetshort+1:end,:);orderedset(1:offsetshort,:)];
shortdiff=orderedsetoffsetshort-orderedset;
shortgrad=atan2(shortdiff(:,2),shortdiff(:,1));   %angle in radians
shortgradoffset=[shortgrad(gradientoffsetshort+1:end,:);shortgrad(1:gradientoffsetshort,:)];
shortgraddiff=shortgradoffset-shortgrad;
shortgraddiff=shortgraddiff+2*pi*(shortgraddiff<0);  %account for 4 quadrants
%%%%%%% long steps %%%%%%%%%%%%%%%%%%%
orderedsetoffsetlong=[orderedset(offsetlong+1:end,:);orderedset(1:offsetlong,:)];
longdiff=orderedsetoffsetlong-orderedset;
longgrad=atan2(longdiff(:,2),longdiff(:,1));
longgradoffset=[longgrad(gradientoffsetlong+1:end,:);longgrad(1:gradientoffsetlong,:)];
longgraddiff=longgradoffset-longgrad;
longgraddiff=longgraddiff+2*pi*(longgraddiff<0);

%%% find deflections above threshold %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% short steps %%%%%%%%%%%%%%%%%%
shortgraddiff(shortgraddiff>=pi)=0;  %exclude
vIdxmaskshort=shortgraddiff>shortgradthresh;
%%%%%%% long steps %%%%%%%%%%%%%%%%%%%
vIdxmasklong=longgraddiff>longgradthresh & longgraddiff<pi;
vIdxmasklong=[zeros(offsetlong,1);vIdxmasklong(1:end-offsetlong)];
vIdxmasklong=imdilate(vIdxmasklong,strel('square',1+nucr/2));

%%% find local maxima of short steps %%%%%%%%%%%%%%%%%%%%%%%
vIdxmaskshort=imclose(vIdxmaskshort,strel('square',3));  %join proximal deflection islands
vIdxobs=regionprops(bwlabel(vIdxmaskshort),'PixelIdxList');
maxmask=zeros(size(vIdxmaskshort));
for rpc=1:length(vIdxobs)
    pix=vIdxobs(rpc).PixelIdxList;
    [~,index]=max(shortgraddiff(pix));
    maxmask(pix(index)+offsetshort)=1;
end
maxmask=maxmask(1:numpix);  %remove any overhang

%%% find coincidence of long mask & local maxima of short mask %%%%%%%%%%%%
vIdxmask=vIdxmasklong & maxmask;
vIdx=find(vIdxmask);

%%% mark candidate vertices [debugging purposes] %%%%%%%%%%%%%%%%%%%%%%%%%%
if bridgeimage
    vpos=orderedset(vIdx,:);
    for vc=1:size(vpos,1)
        bordermask(vpos(vc,2),vpos(vc,1))=1;
    end
end

%%% pair and connect vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vnum=length(vIdx);
if vnum>=2
    periIdx=vIdx;
    perisize=length(vIdxmask);
    periIdxadj1=[periIdx(2:end);perisize+periIdx(1)];
    pairperi1=periIdxadj1-periIdx;    %perimeter distance between adj vertices
end
while vnum>=2
    skipvertices=0;
    vpos=orderedset(vIdx,:);
    vposadj1=[vpos(2:end,:);vpos(1,:)];
    pair1=vposadj1-vpos;
    pairdist1=sqrt(sum(pair1.^2,2));
    [bestpair1,ordx1]=sort(pairperi1./pairdist1);
    idx=ordx1(end);
    if idx==vnum
        idxadj=1;
    else
        idxadj=idx+1;
    end
    if vnum>=5
        vposadj2=[vposadj1(2:end,:);vposadj1(1,:)];
        pair2=vposadj2-vpos;
        pairdist2=sqrt(sum(pair2.^2,2));
        pairperi2=[pairperi1(2:end);pairperi1(1)];
        pairperi2=pairperi1+pairperi2;    %perimeter btwn every other vertice
        [bestpair2,ordx2]=sort(pairperi2./pairdist2);
        if bestpair2(end)>bestpair1(end)
            skipvertices=1;
            idx=ordx2(end);
            if idx==vnum
                idxadj=2; idxint=1;
            elseif idx==vnum-1
                idxadj=1; idxint=vnum;
            else
                idxadj=idx+2; idxint=idx+1;
            end
        end
    end

    if vnum<=3
        if vnum==3
            idxadjadj=idxadj+1;
            if idxadjadj>vnum
                idxadjadj=1;
            end
            pairperi1(idxadj)=pairperi1(idxadj)+pairperi1(idxadjadj);
            pairperi1(idxadjadj)=nucr/2*pi;
        end
        vlo=pairperi1<nucr/2*pi;   %less than half circumference of avg size nucleus
        if sum(vlo)             %either side too short
            break
        end
    end

    [bcx,bcy]=bridge(vpos(idx,:),vpos(idxadj,:));
    for bci=1:length(bcx)
        bridgeflag=1;
        bordermask(bcy(bci),bcx(bci))=1;
    end
    
    %%% assign new perimter distance & remove old vertices %%%%%%%%%%%%%%%%
    previdx=idx-1;
    if previdx==0
        previdx=vnum;
    end
    pairperi1(previdx)=pairperi1(previdx)+length(bcx)-1+pairperi1(idxadj);
    if skipvertices==0
        vIdx([idx,idxadj])=[];
        pairperi1([idx,idxadj])=[];
    else
        vIdx([idx,idxint,idxadj])=[];
        pairperi1([idx,idxint,idxadj])=[];
    end
    vnum=length(vIdx);
end
end