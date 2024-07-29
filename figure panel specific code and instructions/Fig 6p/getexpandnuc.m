function exp_nuc=getexpandnuc(nuc_label,ringwidth,raw)

% Uses a mask and ringwidth as an input to creat

%%% define cytosolic ring radii %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(nuc_label);
nuc_mask=nuc_label>0;
%%% define border margins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ringzeroborder=ones(height,width);
ringzeroborder(1:ringwidth+1,:)=0;ringzeroborder(end-ringwidth:end,:)=0;ringzeroborder(:,1:ringwidth+1)=0;ringzeroborder(:,end-ringwidth:end)=0;

%%% define midband of cytoring with same labels as nuclei %%%%%%%%%%%%%%%%%
ringedgeouter=imdilate(nuc_label,strel('disk',ringwidth,0));  %outerrad one pixel greater than midrad
%ringedgeinner=imdilate(nuc_label,strel('disk',midrad-offset,0));  %outerrad one pixel greater than midrad
%cytoring=ringedgeouter-nuc_label;                     %define cytoring with label ID
%cytoring=ringedgeouter-ringedgeinner;

%%% detect cell-cell borders and inflate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
borders=bwmorph(nuc_label,'bothat');
% imdilated=ringedgeouter>0;
% imclosed=imerode(imdilated,strel('disk',ringwidth,0));
% borders=imclosed-nuc_mask;
borders=imdilate(borders,strel('disk',2,0));

%%% clarify boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_nuc=ringedgeouter.*~borders;
exp_nuc=exp_nuc.*ringzeroborder;

%%%%%%%% reconstitute absent rings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cru=unique(ringedgeouter);
fcru=unique(exp_nuc);
noring=cru(~ismember(cru,fcru));
cytoring_rzb=ringedgeouter.*ringzeroborder;
if ~isempty(noring)
    for i=noring'
        %ring_label(cytoring_rzb==i)=i;
        exp_nuc(ringedgeouter==i)=i;
    end
end



%% visualization for debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tempframe=imadjust(mat2gray(raw));
% tempframe(:,:,2)=bwmorph(ringedgeouter,'remove');
% tempframe(:,:,3)=bwmorph(nuc_mask,'remove');
% tempframe(:,:,3)=bwmorph(borders,'remove');
% tempframe(:,:,2)=exp_nuc>0;
% tempframe(:,:,3)=0;
% imshow(tempframe);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end