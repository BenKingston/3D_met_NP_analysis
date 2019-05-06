


function [pre_nuclei, pre_vessels, pre_micromet] = pre_process_func(nuclei,vessels,micromets,save_dir,sample_name)

shortfile = sample_name;
display(['Pre-processing ' shortfile])

%Pre-processing of nuclei channel
nuclei_single = single(nuclei);
nuclei_local = nuclei_single./(0.2*max(nuclei_single(:))+(im2mat(gaussf(nuclei_single,size(nuclei_single,1)/10))));
nuclei_loglocal = mat2im(log(nuclei_local+0.1));

nuclei_processed = uint8(stretch(nuclei_loglocal));
pre_nuclei = nuclei_processed;

%Pre-processing of blood vessels channel
vessels_single = single(vessels);
vessels_local = vessels_single./(0.2*max(vessels_single(:))+(im2mat(gaussf(vessels_single,size(vessels_single,1)/10))));
vessels_loglocal = mat2im(log(vessels_local+0.1));

vessels_processed = uint8(stretch(vessels_loglocal));
pre_vessels = vessels_processed;

%Pre-processing of ki67 channel
ki67_single = single(micromets);
ki67_local = ki67_single./(0.2*max(ki67_single(:))+(im2mat(gaussf(ki67_single,size(ki67_single,1)/10))));
ki67_loglocal = mat2im(log(ki67_local+0.1));

ki67_processed = uint8(stretch(ki67_loglocal));
pre_micromet = ki67_processed;

%Write pre-processed files
cd(save_dir)

nuclei_processed_name = strcat(shortfile,'_pre_processed_nuclei.tif');
vessels_processed_name = strcat(shortfile,'_pre_processed_vessels.tif');
ki67_processed_name = strcat(shortfile,'_pre_processed_micromets.tif');

num_slices = size(nuclei,3);

imwrite(uint8(nuclei_processed(:,:,1)),nuclei_processed_name);
        
   for p = 2:num_slices
            imwrite(uint8(nuclei_processed(:,:,p)),nuclei_processed_name, 'WriteMode','append');
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
