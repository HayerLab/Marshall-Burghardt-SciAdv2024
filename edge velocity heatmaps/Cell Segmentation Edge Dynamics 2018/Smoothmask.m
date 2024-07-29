close all;
%mask_raw=MaskThisCell{51};

for frameNum=1:90
    mask_raw=MaskThisCell{frameNum};
    se=strel('disk',4);
    mask_shrink=imerode(mask_raw,se);
    mask_dil=imdilate(mask_shrink,se);
    BWoutline=bwperim(mask_dil);
    
    imtemp=imFRETOutline{frameNum};
    redchan=imtemp(:,:,1);
    greenchan=imtemp(:,:,2);
    bluechan=imtemp(:,:,3);
    redchan(BWoutline)=1;
    greenchan(BWoutline)=0;
    bluechan(BWoutline)=0;
    imtemp(:,:,1)=redchan;
    imtemp(:,:,2)=greenchan;
    imtemp(:,:,3)=bluechan;

    imwrite(imtemp,'test.tiff','WriteMode','append','Compression','none');
end
%%