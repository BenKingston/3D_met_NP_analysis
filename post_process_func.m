


function [post_nuclei, post_nuclei_dilate, post_vessels, post_micromet] = post_process_func(pre_nuclei,seg_nuclei,seg_vessels,seg_micromets,save_dir,sample_name)

shortfile = sample_name;
display(['Post-processing ' shortfile])

    dapi = pre_nuclei;

    fill_tissue = fillholes(dapi>40);
    erode_tissue = erosion(fill_tissue,14,'elliptic');
    open_tissue = opening(erode_tissue,10,'elliptic');
    gauss_tissue = gaussf(open_tissue,3);
    block_tissue = label(gauss_tissue>0.5,1,50000000,0); 
    tissue_area = uint8(block_tissue>0);
    
    num_slices = size(tissue_area,3);

%%%%Load in Ilastik segmented image files
nuclei_seg = seg_nuclei;
vessels_seg = seg_vessels;
ki67_seg = seg_micromets;

%%%%%Post-processing of Ilastik segmented nuclei channel
nuclei_seg_bin = nuclei_seg==1;
tissue_area_bin = tissue_area>0;
nuclei_seg_crop = (nuclei_seg_bin.*tissue_area_bin)>0;


threshnuc = nuclei_seg_crop>0;
    threshnuc3 = opening(threshnuc,4);
    dt_threshnuc = dt(threshnuc3);
    
    seeds = maxima(dt_threshnuc,2,0);
    seeds2 = dilation(seeds,4.5)>0;
    image_out = waterseed(seeds2,max(dt_threshnuc)-dt_threshnuc,1,0,0);
    
    threshnuc4 = (uint8(threshnuc3))>0;
    threshnuc4_dil_se = strel('sphere',1);
    threshnuc4_dil = imdilate(threshnuc4,threshnuc4_dil_se);
    threshnuc4_dil(image_out) = false; 
    
nuclabel = label(threshnuc4_dil>0,1,25,0);
nuc_dilated_label = dip_growregions(nuclabel,[],[],1,2,'low_first');
final_nuc = nuclabel;
final_dilated_nuc = nuc_dilated_label;

nuclei_processed = uint32(final_nuc);
nuclei_dilated_processed = uint32(final_dilated_nuc);
%%%%%Post-processing of Ilastik segmented blood vessel channel
vessel_seg_bin = vessels_seg==1;
vessel_seg_crop = vessel_seg_bin.*tissue_area_bin;
se_erode = strel('sphere',2);
vessel_seg_crop_er = imerode(vessel_seg_crop,se_erode);

vessels_processed = vessel_seg_crop_er>0;

%%%%%Post-processing of Ilastik segmented ki67/metastases channel
ki67_bin = ki67_seg==1;
ki67_seg_crop = ki67_bin.*tissue_area_bin;
se_open_ki67 = strel('sphere',3);
ki67_open = imopen(ki67_seg_crop,se_open_ki67);

ki67_label = label(ki67_open>0,1,30000,0);
ki67_processed = uint16(ki67_label);

%%%%%Final segmentation images
post_nuclei = nuclei_processed;
post_vessels = vessels_processed;
post_micromet = ki67_processed;
post_nuclei_dilate = nuclei_dilated_processed;
%%%%%Write post-processed files

nuclei_processed_name = strcat(shortfile,'_post_processed_nuclei.tif');
nuclei_dilated_processed_name = strcat(shortfile,'_post_processed_dialted_nuclei.tif');
vessels_processed_name = strcat(shortfile,'_post_processed_vessels.tif');
ki67_processed_name = strcat(shortfile,'_post_processed_ki67.tif');

num_slices = size(vessels_seg,3);

cd(save_dir)
imwrite(uint16(nuclei_processed(:,:,1)),nuclei_processed_name);
         
  for p = 2:num_slices
          imwrite(uint16(nuclei_processed(:,:,p)),nuclei_processed_name, 'WriteMode','append');
  end

imwrite(uint16(nuclei_dilated_processed(:,:,1)),nuclei_dilated_processed_name);
         
  for p = 2:num_slices
          imwrite(uint16(nuclei_dilated_processed(:,:,p)),nuclei_dilated_processed_name, 'WriteMode','append');
  end  

 imwrite(uint8(vessels_processed(:,:,1)),vessels_processed_name);
         
    for p = 2:num_slices
            imwrite(uint8(vessels_processed(:,:,p)),vessels_processed_name, 'WriteMode','append');
    end   

 imwrite(uint8(ki67_processed(:,:,1)),ki67_processed_name);
         
    for p = 2:num_slices
             imwrite(uint8(ki67_processed(:,:,p)),ki67_processed_name, 'WriteMode','append');
    end   

end


