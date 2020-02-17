function [varargout]=MamutXML2Matlab(varargin)
%% explanation for the input and output
%the input should be 1= the xml file 2= time between each frame
%the output is in order a 1.full list of spots and what track they belong to,
%2. the full spot list, 3. the track_list


%% read the xml file and extract the spots and tracks info
%creates a progress bar
f = uifigure;
if length(varargin) > 1
    app_pos = varargin{2};
   f.Position = [app_pos(1)+0.5*app_pos(3)-200 app_pos(2)+0.5*app_pos(4)-50 400 100];
end

d = uiprogressdlg(f,'Title','Please wait','Message','Reading the xml-file');
d.Value = 0.1;

%read the xml-file
XMLstruct=xml2struct(varargin{1});
all_spots=XMLstruct.TrackMate.Model.AllSpots;
all_tracks=XMLstruct.TrackMate.Model.AllTracks;

%first set the default scaling
x_scale = 1;
y_scale = 1;
z_scale = 1;
scale_unit = XMLstruct.TrackMate.Model.Attributes.spatialunits;
%investigate whether there is a calibration in xml file
if XMLstruct.TrackMate.Model.Attributes.spatialunits == "pixels"
    promt = {'Size of X:','Size of Y:','Size of Z','Unit:'};
    dlgtitle = 'Calibrate scaling';
    dims = [1 30];
    definpu = {'1','1','1','pixels'};
    scaling = inputdlg(promt,dlgtitle,dims,definpu);
    if isempty(scaling)
        scaling = definpu;
    end
    while isnan(str2double(scaling{1})) || isnan(str2double(scaling{2})) || isnan(str2double(scaling{3}))
        promt = {'Size of X (most be number):','Size of Y (most be a number):','Size of Z','Unit:'};
        dlgtitle = 'Calibrate scaling, cancling will set scale to 1 pixel';
        scaling = inputdlg(promt,dlgtitle,dims,scaling);
        if isempty(scaling)
            scaling = definpu;
        end
    end
    x_scale = str2double(scaling{1});
    y_scale = str2double(scaling{2});
    z_scale = str2double(scaling{3});
    scale_unit = scaling{4};
end


d.Value = 0.4;
d.Message = 'Converting to table';

%deal with the spots
num_row=str2double(all_spots.Attributes.nspots);
num_col=9;
col_names=["name","ID", "pos_x","pos_y","pos_z","pos_T","frame","radius","Quality"];
col_types=["string","double","double","double","double","double","double","double","double"];
spots_table=table('size',[num_row,num_col],'VariableNames', col_names, 'variableTypes', col_types);
counts=1;
%loop over the number of Frames
for k=1:length(all_spots.SpotsInFrame)
    %loop over the spots in each frame
    for i=1:length(all_spots.SpotsInFrame{k}.Spot)
        %assing values to the table
        spots_table{counts,1:9}=[string(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.name),...
                                str2double(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.ID),...
                                str2double(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.POSITION_X)*x_scale,...
                                str2double(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.POSITION_Y)*y_scale,...
                                str2double(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.POSITION_Z)*z_scale,...
                                str2double(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.POSITION_T),...
                                str2double(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.FRAME),...
                                str2double(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.RADIUS),...
                                str2double(all_spots.SpotsInFrame{k}.Spot{i}.Attributes.QUALITY)];
        counts=counts+1;
    end
end

%deal with the tracks
%determine the number of rows needed:
num_row=0;
for k=1:length(all_tracks.Track)
    num_row=num_row+length(all_tracks.Track{k}.Edge);
end
num_col=5;
col_names=["Track_ID","Spot_source_ID","Spot_target_ID","Track_name","Displacement"];
col_types=["double","double","double","string","double"];
track_table=table('size',[num_row,num_col],'VariableNames', col_names, 'variableTypes', col_types);
counts=1;
%assign value to the table
for k=1:length(all_tracks.Track)
    track_table{counts:counts+(length(all_tracks.Track{k}.Edge)-1),[1,4]}=repmat([str2double(all_tracks.Track{k}.Attributes.TRACK_ID),string(all_tracks.Track{k}.Attributes.name)],length(all_tracks.Track{k}.Edge),1);
    for i=1:length(all_tracks.Track{k}.Edge)
        track_table{counts,[2,3,5]}=[str2double(all_tracks.Track{k}.Edge{i}.Attributes.SPOT_SOURCE_ID),...
                                    str2double(all_tracks.Track{k}.Edge{i}.Attributes.SPOT_TARGET_ID),...
                                    str2double(all_tracks.Track{k}.Edge{i}.Attributes.DISPLACEMENT)];
        counts=counts+1;
    end
end


d.Value = 0.6;
d.Message = 'Merging spots and tracks';

%create a table that cotain the spots but also the which track the spots
%belong to.
spots_tracks = spots_table;
spots_tracks.TrackID = zeros(height(spots_table),1);
spots_tracks.SourceID = zeros(height(spots_table),1);


for i=1:height(spots_tracks)
    %the max is there to ensure that only one value is transfered (wheresas
    %any spot with both a source and target would give 2 identical numbers)
    track_ID = max(table2array(track_table(track_table.Spot_source_ID == spots_tracks.ID(i) | track_table.Spot_target_ID == spots_tracks.ID(i),1)));
    sourceID = track_table.Spot_source_ID(track_table.Spot_target_ID == spots_tracks.ID(i));
    targetID = track_table.Spot_target_ID(track_table.Spot_source_ID == spots_tracks.ID(i));

    if ~isempty(track_ID)
        spots_tracks.TrackID(i)=track_ID;
    else
        spots_tracks.TrackID(i)=-1;
    end
    
    if ~isempty(sourceID)
        if length(sourceID) > 1
            sourceID = max(sourceID);
        end
        spots_tracks.SourceID(i) = sourceID;
    else
        spots_tracks.SourceID(i) = -1;
    end
    
    if ~isempty(targetID)
        if length(targetID) > 1
            targetID = max(targetID);
        end
        spots_tracks.TargetID(i) = targetID;
    else
        spots_tracks.TargetID(i) = -1;
    end
    
end
%% Deal with daughter cells to determine which sister cells belong to which
%sister cell

d.Value = 0.8;
d.Message = 'Getting siblings to talk to each other';

spots_tracks.TrackID_sub = string(blanks(height(spots_table))');
uni_tracks=unique(spots_tracks.TrackID);

for i=1:length(uni_tracks)
    info_temp = spots_tracks(uni_tracks(i) == spots_tracks.TrackID,:);
    %test whether ther is any cell division
    if height(info_temp) > length(unique(info_temp.frame))
        [~,I]=sort(info_temp.frame);
        info_temp=info_temp(I,:);
        info_temp.TrackID_sub(1) = string(num2str(info_temp.TrackID(1)));

        %when cell division is detected, determine where they are and make
        %sub name for them, after a division they will get the name of the
        %mother + A & B
        for k=1:(height(info_temp)-1)
            Targets = info_temp(info_temp.SourceID == info_temp.ID(k),:);
            Trackname = info_temp.TrackID_sub(k);
            spots_tracks.TrackID_sub(spots_tracks.ID == info_temp.ID(k)) = info_temp.TrackID_sub(k);
            
            if height(Targets) == 2
                info_temp.TrackID_sub(info_temp.SourceID == info_temp.ID(k)) = string([char(Trackname) 'A';char(Trackname) 'B']);
            elseif height(Targets) > 2
                disp("multi split!?")
            else
               info_temp.TrackID_sub(info_temp.SourceID == info_temp.ID(k)) =  info_temp.TrackID_sub(k);
            end
                       
        end
        spots_tracks.TrackID_sub(spots_tracks.ID == info_temp.ID(height(info_temp))) = info_temp.TrackID_sub(height(info_temp));
    else
        spots_tracks.TrackID_sub(spots_tracks.TrackID == uni_tracks(i)) = string(num2str(uni_tracks(i)));
    end
    
    
end

d.Value = 1;
d.Message = 'Performing house keeping';
pause(1)
close(d)
close(f)

%% Output (temp for now)
varargout{1}=spots_tracks;
varargout{2}=spots_table;
varargout{3}=track_table;
varargout{4}={x_scale,y_scale,z_scale,scale_unit};
