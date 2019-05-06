
function [int_file] = convert_czi_to_ometiff_revised_func(data_folder)

matlab_folder = pwd;%remembers the current folder so it can come back to it
%load required variables

%data_folder = uigetdir('','Select folder with image files');

cd(data_folder); 
% answer = inputdlg('How many hours from now should files be converted?');
% disp('delayed for')
% disp(strcat(cell2mat(answer), ' hours'))
% 
% pause(str2double(cell2mat(answer))*60*60);

files = dir('*.czi');    

%%
for int_file = 1:size(files,1)   
%%
    [~,shortfile] = fileparts(files(int_file).name);
    newimgdir = strcat(data_folder,'\',shortfile);
    mkdir(newimgdir);                            
    cd(newimgdir);                               
    
    origfilename_ch1 = strcat(shortfile, '_iso_ch1.tif'); 
    origfilename_ch2 = strcat(shortfile, '_iso_ch2.tif');
    origfilename_ch3 = strcat(shortfile, '_iso_ch3.tif');
    origfilename_ch4 = strcat(shortfile, '_iso_ch4.tif');
    origfilename_ch5 = strcat(shortfile, '_iso_ch5.tif');
    filenames = [origfilename_ch1; origfilename_ch2; origfilename_ch3; origfilename_ch4; origfilename_ch5];  
    
    metadatafilename = strcat(shortfile, '_iso_info.csv');
    metadatafullfilename = strcat(shortfile, '_iso_fullmeta.csv'); 
    %%
    if exist(metadatafullfilename,'file')==2    
        disp('skipping');
        disp(origfilename_ch1);
    else
        
        
        %tic
        cd(data_folder);
			reader = bfGetReader(files(int_file).name);
			omeMeta = reader.getMetadataStore();  
        % data = bfopen(files(int_file).name);       
        %disp('reading time');
        %toc
%         tic
       cd(newimgdir);
    
        
        stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
        stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
        stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
        numchannels = omeMeta.getPixelsSizeC(0).getValue();
        timepoints = omeMeta.getPixelsSizeT(0).getValue();
        
        metadata = reader.getGlobalMetadata();
        metadataKeys = metadata.keySet().iterator();
        metadataarray = {};
        for m=1:metadata.size()
            key = metadataKeys.nextElement();
            value = metadata.get(key);
            %   fprintf('%s = %s\n', key, value)
            metadataarray = [metadataarray; {key} {value}];
        end
       
        [metadatakeys_cell, I]= sort( metadataarray(:,1));
        metadatavalues_cell =  metadataarray(I,2);
        
        %%
        T2 = table(metadatakeys_cell,metadatavalues_cell);
        [laser_wavelength, laser_powers, cam_int_time, filterset, ...
            trackconfig, leftrightalignment, LSthickness, imagezoom, ... 
            objective, views] = extractMetaData(metadatakeys_cell,metadatavalues_cell);
        %         strphysSizeX = (metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'));
        %         strphysSizeY = (metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1'));
        %         strphysSizeZ = (metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1'));
        
        %%
        physSizeX = str2num(metadata.get('Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'));
        physSizeY = str2num(metadata.get('Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1'));
        physSizeZ = str2num(metadata.get('Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1'));
        
        newSizeX = round(stackSizeX*physSizeX/physSizeZ);
        newSizeY = round(stackSizeY*physSizeY/physSizeZ);
        newSizeZ = stackSizeZ;
        
        newphysX = physSizeX*stackSizeX/newSizeX;
        newphysY = physSizeY*stackSizeY/newSizeY;
        newphysZ = physSizeZ;
        plane = zeros(newSizeX,newSizeY,newSizeZ,numchannels,'uint16');
        %%
  
        %%
        
         for k = 1:numchannels
            for p = 1:newSizeZ
                series1_plane1 = bfGetPlane(reader, p*numchannels + k - numchannels);
                plane(:,:,p,k) = imresize(series1_plane1, [newSizeX newSizeY]);
                waitbar((p+newSizeZ*(k-1))/(numchannels*newSizeZ))
            end
         end

%%        
        names = {'pixel size X'; 'pixel size Y'; 'pixel size Z'; ...
            '# of Pixels X'; '# of Pixels Y'; '# of Pixels Z'; ...
            'old pixel size'; 'Total Channels'; 'Vessel Channel'; ...
            'Particle Channel'; 'Zoom'; 'LightSheet Thickness'; ...
            'Objective'; 'Laser 1 nm'; 'Laser 2 nm'; 'Laser 3 nm';...
            'Laser 4 nm'; 'Laser 5 nm'; 'Laser Power 1';...
            'Laser Power 2'; 'Laser Power 3';...
            'Laser Power 4'; 'Laser Power 5'; 'Filters 1'; 'Filters 2';...
            'Filters 3'; 'Filters 4'; 'Filters 5'; 'Exposure 1'; ...
            'Exposure 2'; 'Exposure 3'; 'Exposure 4'; 'Exposure 5';...
            'Track1'; 'Track2'; 'Track3'; 'Track4'; 'Track5'; 'ZLeft.405';...
            'ZLeft.488'; 'ZLeft.561'; 'ZLeft.638'; 'ZRight.405'; ...
            'ZRight.488'; 'ZRight.561'; 'ZRight.638' };
        newphys = {newphysX; newphysY; newphysZ; newSizeX; newSizeY; ...
            newSizeZ; physSizeX; numchannels; 3; numchannels; ...
            imagezoom; LSthickness; objective{1}; laser_wavelength(1);...
            laser_wavelength(2); laser_wavelength(3); laser_wavelength(4);...
            laser_wavelength(5); laser_powers(1); laser_powers(2); ...
            laser_powers(3); laser_powers(4); laser_powers(5);...
            filterset{1}; filterset{2}; filterset{3};...
            filterset{4}; filterset{5}; cam_int_time(1); cam_int_time(2);...
            cam_int_time(3);cam_int_time(4);cam_int_time(5); ...
            trackconfig{1};trackconfig{2};trackconfig{3};trackconfig{4};...
            trackconfig{5};leftrightalignment(1);leftrightalignment(2);...
            leftrightalignment(3);leftrightalignment(4);...
            leftrightalignment(5);leftrightalignment(6);...
            leftrightalignment(7);leftrightalignment(8);};
        
        T = table(names, newphys );
        cd(newimgdir);
        writetable(T,metadatafilename);
        writetable(T2,metadatafullfilename);
        disp('processing time');
%         toc
        
                tic
        for k = 1:numchannels
           cd(newimgdir);
            clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(plane(:,:,:,k), filenames(k,:), options);

        end
        reader.close()
        clear reader
        disp('writing time')
        toc
     

        
        %         save(metadatafilename,'stackSizeX','stackSizeY','stackSizeZ','numchannels','physSizeX','physSizeY','physSizeZ','newSizeX','newSizeY','newSizeZ','newphysX','newphysY','newphysZ')
    end
    cd(data_folder);
    
    
    delete(files(int_file).name);
end

cd(matlab_folder)



end