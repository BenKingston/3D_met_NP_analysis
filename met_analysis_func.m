

function [met_cell_np_int, met_cell_dist, met_dilate_nuclei_label] = met_analysis_func(pre_nuclei,post_nuclei,post_vessels,post_micromets,nanoparticle_ch,save_dir,sample_name)


shortfile = sample_name;
display(['Processing ' shortfile])
tic

dapi = pre_nuclei;

    fill_tissue = fillholes(dapi>40);
    erode_tissue = erosion(fill_tissue,14,'elliptic');
    open_tissue = opening(erode_tissue,10,'elliptic');
    gauss_tissue = gaussf(open_tissue,3);
    block_tissue = label(gauss_tissue>0.5,1,50000000,0); 
    tissue_area = uint8(block_tissue>0);


nuclei_post_processed = post_nuclei;
vessels_post_processed = post_vessels;
ki67_post_processed = post_micromets;

%%% Load particle channel and trim to tissue boundary
particle = nanoparticle_ch;

tissue_area_bin = tissue_area>0;
particle_crop = particle.*uint16(tissue_area_bin);
ki67_notmet_bin = ki67_post_processed <1;

background_mask = particle_crop > 3000;
background_img = uint16(background_mask).*particle_crop;
background_img2 = uint16(ki67_notmet_bin).*background_img;

num_slices = size(ki67_post_processed,3);

%%% Analysis of individual mets
nuclabel = label(nuclei_post_processed>0,1,25,0);
nuc_dilated_label = dip_growregions(nuclabel,[],[],1,2,'low_first');
nuc_dilated_label = uint32(nuc_dilated_label);

Ind_met_stats = regionprops3(ki67_post_processed,particle_crop,'Centroid','MeanIntensity','Volume','SurfaceArea');

ki67_post_processed_bin = ki67_post_processed>0;

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

    
new_all_nuc_int_crop = new_all_nuc_int.*uint32(ki67_post_processed_bin);    
new_all_nuc_dist_crop = new_all_nuc_dist.*uint32(ki67_post_processed_bin);      
nuclei_dilated_post_processed_crop = nuc_dilated_label.*uint32(ki67_post_processed_bin);      

%%% Single cell analysis of individual mets

max_metnuc = max(ki67_post_processed(:));
min_metnuc = min(ki67_post_processed(:))+1;
    for a = min_metnuc:max_metnuc
        show_met = ki67_post_processed;
        show_met(show_met~=a) = 0;
        isolate_met = (show_met>0);
        nucinmet = uint32(isolate_met).*nuclei_dilated_post_processed_crop;
        

        msr_met_cell_particle_intensity = regionprops3(uint32(nucinmet),particle_crop,'MeanIntensity','Centroid');
        msr_met_cell_dist_vess = regionprops3(uint32(nucinmet),dt_vessels,'MeanIntensity','Centroid');
        
        msr_met_cell_particle_intensity_nonan = rmmissing(msr_met_cell_particle_intensity);
        msr_met_cell_dist_vess_nonan = rmmissing(msr_met_cell_dist_vess);
        
        msr_met_cell_particle_intensity_nonan.Properties.VariableNames = {'CentroidNPIntensity' 'MeanNPIntensity'};
        msr_met_cell_dist_vess_nonan.Properties.VariableNames = {'CentroidDistIntensity' 'MeanDistIntensity'};
        
        
        Met_cell_npint_dist = [msr_met_cell_particle_intensity_nonan msr_met_cell_dist_vess_nonan];
        
        Met_t_name_nuc = strcat(shortfile,'-Metstat-Met-_np_',num2str(a),'.csv');
    cd(save_dir)
    writetable(Met_cell_npint_dist,Met_t_name_nuc);
    
    end

met_cell_np_int = new_all_nuc_int_crop;
met_cell_dist = new_all_nuc_dist_crop;
met_dilate_nuclei_label = nuclei_dilated_post_processed_crop;
    
%%% Writing tables and images
cd(save_dir)

all_met_table_name = strcat(shortfile,'all-Metstat_np','.csv');    
writetable(Ind_met_stats, all_met_table_name);

met_mean_np_int_name = strcat(shortfile,'Met mean NP Intensity','.tif');
met_mean_dist_name = strcat(shortfile,'Met mean distance to vessel','.tif');
met_nuclei_label_name = strcat(shortfile,'Met labeled cells','.tif');

num_slices = size(ki67_post_processed,3);

imwrite(uint16(new_all_nuc_int_crop(:,:,1)),met_mean_np_int_name);
         
  for p = 2:num_slices
          imwrite(uint16(new_all_nuc_int_crop(:,:,p)),met_mean_np_int_name, 'WriteMode','append');
  end

imwrite(uint16(new_all_nuc_dist_crop(:,:,1)),met_mean_dist_name);
         
  for p = 2:num_slices
          imwrite(uint16(new_all_nuc_dist_crop(:,:,p)),met_mean_dist_name, 'WriteMode','append');
  end

imwrite(uint16(nuclei_dilated_post_processed_crop(:,:,1)),met_nuclei_label_name);
         
  for p = 2:num_slices
          imwrite(uint16(nuclei_dilated_post_processed_crop(:,:,p)),met_nuclei_label_name, 'WriteMode','append');
  end
toc
end