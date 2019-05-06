function [met_cell_np_int, met_cell_dist, met_dilate_nuclei_label] = tumour_analysis_func(pre_nuclei,post_nuclei,post_vessels,nanoparticle_ch,save_dir,sample_name)


shortfile = sample_name;
display(['Processing ' shortfile])
tic

%Trim the tissue boundary using the nuclei channel
dapi = pre_nuclei;

    fill_tissue = fillholes(dapi>40);
    erode_tissue = erosion(fill_tissue,14,'elliptic');
    open_tissue = opening(erode_tissue,10,'elliptic');
    gauss_tissue = gaussf(open_tissue,3);
    block_tissue = label(gauss_tissue>0.5,1,50000000,0); 
    tissue_area = uint8(block_tissue>0);


nuclei_post_processed = post_nuclei;
vessels_post_processed = post_vessels;

%%% Load particle channel and trim to tissue boundary

particle = nanoparticle_ch;

tissue_area_bin = tissue_area>0;
particle_crop = particle.*uint16(tissue_area_bin);

%%% Analysis of tumour tissue
nuclabel = label(nuclei_post_processed>0,1,25,0);
nuc_dilated_label = dip_growregions(nuclabel,[],[],1,2,'low_first');
nuc_dilated_label = uint32(nuc_dilated_label);

Ind_tumor_stats = regionprops3(tissue_area_bin,particle_crop,'Centroid','MeanIntensity','Volume','SurfaceArea');

all_cells_int_stats = regionprops3(nuc_dilated_label,particle_crop,'MeanIntensity','VoxelIdxList');
all_cells_int_stats_struct = table2struct(all_cells_int_stats);

vess_thresh = vessels_post_processed>0;
dt_vessels = bwdist(vess_thresh);
all_cells_dist_stats = regionprops3(nuc_dilated_label,dt_vessels,'MeanIntensity','VoxelIdxList');
all_cells_dist_stats_struct = table2struct(all_cells_dist_stats);

max_cell_in_tissue = max(nuc_dilated_label(:));
min_cell_in_tissue = min(nuc_dilated_label(:))+1;

new_all_nuc_int =  uint32(nuc_dilated_label);     
    
    for b = min_cell_in_tissue:max_cell_in_tissue
            new_all_nuc_int(all_cells_int_stats_struct(b).VoxelIdxList) = all_cells_int_stats_struct(b).MeanIntensity;
    end

new_all_nuc_dist =  uint32(nuc_dilated_label);     

    for b = min_cell_in_tissue:max_cell_in_tissue
            new_all_nuc_dist(all_cells_dist_stats_struct(b).VoxelIdxList) = all_cells_dist_stats_struct(b).MeanIntensity;
    end
    
new_all_nuc_int_crop = new_all_nuc_int;    
new_all_nuc_dist_crop = new_all_nuc_dist;      
nuclei_dilated_post_processed_crop = nuc_dilated_label;      

met_cell_np_int = new_all_nuc_int_crop;
met_cell_dist = new_all_nuc_dist_crop;
met_dilate_nuclei_label = nuclei_dilated_post_processed_crop;


%%% Writing tables and images
cd(save_dir)
all_cells_int_stats.Properties.VariableNames = {'VoxelIdxList_NPint' 'MeanIntensity_NPint'};
all_cells_dist_stats.Properties.VariableNames = {'VoxelIdxList_dist' 'MeanIntensity_dist'};

Tumor_cell_npint_dist = [all_cells_int_stats.MeanIntensity_NPint all_cells_dist_stats.MeanIntensity_dist];
Tumor_cell_table = array2table(Tumor_cell_npint_dist);
Tumor_t_name_nuc = strcat(shortfile,'-Tumor-cells-all','.csv');
writetable(Tumor_cell_table,Tumor_t_name_nuc);

all_met_table_name = strcat(shortfile,'-Tumor-area-','.csv');    
writetable(Ind_tumor_stats, all_met_table_name);

tumor_mean_np_int_name = strcat(shortfile,'Tumor mean NP Intensity','.tif');
tumor_mean_dist_name = strcat(shortfile,'Tumor mean distance to vessel','.tif');
tumor_nuclei_label_name = strcat(shortfile,'Tumor labeled cells','.tif');

num_slices = size(tissue_area,3);

imwrite(uint16(new_all_nuc_int_crop(:,:,1)),tumor_mean_np_int_name);
         
  for p = 2:num_slices
          imwrite(uint16(new_all_nuc_int_crop(:,:,p)),tumor_mean_np_int_name, 'WriteMode','append');
  end

imwrite(uint16(new_all_nuc_dist_crop(:,:,1)),tumor_mean_dist_name);
         
  for p = 2:num_slices
          imwrite(uint16(new_all_nuc_dist_crop(:,:,p)),tumor_mean_dist_name, 'WriteMode','append');
  end

imwrite(uint16(nuclei_dilated_post_processed_crop(:,:,1)),tumor_nuclei_label_name);
         
  for p = 2:num_slices
          imwrite(uint16(nuclei_dilated_post_processed_crop(:,:,p)),tumor_nuclei_label_name, 'WriteMode','append');
  end
  toc
end