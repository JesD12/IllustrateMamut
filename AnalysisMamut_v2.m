function [output,delta_times]= AnalysisMamut_v2(Data)
            
    %% setup basic data and calculate sister distance
    Framelist = unique(Data.frame);
    Tracklist = unique(Data.TrackID);

    TrackdivisionLogical = false(length(Tracklist),1);
    for k=1:length(Tracklist)
        TrackdivisionLogical(k) = ~(length(Data.frame(Data.TrackID == Tracklist(k))) == length(unique(Data.frame(Data.TrackID == Tracklist(k)))));
    end

    SisterDistance = zeros(length(Framelist),length(Tracklist(TrackdivisionLogical)));
    Data.SisterDistance = zeros(height(Data),1);
    losttrack = zeros(length(Framelist)*length(Tracklist(TrackdivisionLogical)),1);
    for i = Tracklist(TrackdivisionLogical)'
        Data_sub = Data(Data.TrackID == i,:);
        for t = 1:length(Framelist)
            Data_sub_frame = Data_sub(Data_sub.frame == Framelist(t),:);
            if height(Data_sub_frame) <= 1
                SisterDistance(t,Tracklist(TrackdivisionLogical) == i) = 0;
            elseif height(Data_sub_frame) == 2
                SisterDistance(t, Tracklist(TrackdivisionLogical) == i) = sum(diff(table2array(Data_sub_frame(:,3:5))).^2)^0.5;
                Data.SisterDistance(Data.TrackID == i & Data.frame == Framelist(t)) = repmat(sum(diff(table2array(Data_sub_frame(:,3:5))).^2)^0.5,...
                                                                        size(Data.SisterDistance(Data.TrackID == i & Data.frame == Framelist(t))));
            elseif height(Data_sub_frame)> 2
                disp("more than 2 sisters see losttrack variable for affected tracks")
                losttrack(t+(i-1)*length(Framelist)) = i;
            else
%                 disp('they are lost')
%                 disp(height(Data_sub_frame))
            end
            if ~isempty(unique(losttrack))
                assignin('base','losttrack',unique(losttrack));
            end
       end
    end


%% calculate neighbor distance


    Data.NearestNeighbourDistance = zeros(height(Data),1);
    Data.NearestNeighbourID = zeros(height(Data),1);
    % for now just need for first frame I think:

%             for t = 1:length(Framelist)
    SelectionVector = Data.frame == Framelist(1);
    IDs=Data.TrackID(SelectionVector);
    pos_x=(Data.pos_x(SelectionVector));
    pos_y=(Data.pos_y(SelectionVector));
    pos_z=(Data.pos_z(SelectionVector));

    x_dif = pos_x-pos_x';
    y_dif = pos_y-pos_y';
    z_dif = pos_z-pos_z';

    distances=(x_dif.^2+y_dif.^2+z_dif.^2).^0.5;
    distances_nan = (distances./distances).*distances;
    [distances_min,idx] = min(distances_nan);
    Data.NearestNeighbourDistance(SelectionVector) = distances_min;
    Data.NearestNeighbourID(SelectionVector) = IDs(idx); 
%                
%             end

    Data.FirstNeighbourDistance = zeros(height(Data),1);
    Tracklist_sub = Data.TrackID(Data.frame == Framelist(1));
    NearestNeighbour_first = Data.NearestNeighbourID(Data.frame == Framelist(1));
    for t = 1: length(Framelist)
        for k = 1:length(Tracklist_sub)

            idx = Data.frame == Framelist(t) & Data.TrackID == Tracklist_sub(k);
            idx_neighbour = Data.frame == Framelist(t) & Data.TrackID == NearestNeighbour_first(k);

            if max(idx) && max(idx_neighbour)
                distance_first = (((Data.pos_x(idx)-Data.pos_x(idx_neighbour)).^2+...
                                    (Data.pos_y(idx)-Data.pos_y(idx_neighbour)).^2+...
                                    (Data.pos_z(idx)-Data.pos_z(idx_neighbour)).^2).^0.5);

                if size(distance_first) == size(Data.FirstNeighbourDistance(idx))
                    Data.FirstNeighbourDistance(idx) = distance_first;
                elseif size(distance_first) >= size (Data.FirstNeighbourDistance(idx))
                    %fix problem where the first nearest neighbor has split up
                    Data.FirstNeighbourDistance(idx) = min(distance_first);
                end
            elseif max(idx)
                Data.FirstNeighbourDistance(idx) = 0;
            end
        end
    Tracklist_sub_name = unique(Data.TrackID_sub);

    end

%% plot convex_Volume
    Volumes=zeros(length(Framelist),1);
    for t=1:length(Framelist)
        [~,Volumes(t)] = convhull(table2array(Data(Data.frame == Framelist(t),3:5)));
    end
    figure(1);
    clf;
    plot(Volumes)
    title("Volume of Convex Hull")
    xlabel("Frame")
    ylabel("µm^3")

    
%% plot SisterDistancePlot
    if size(SisterDistance) > 0
        colormaps_SD=distinguishable_colors(size(SisterDistance,2));
        figure(2);
        clf;
        assignin('base','sis',SisterDistance)
        plot(1:(length(nonzeros(SisterDistance(:,1)))),nonzeros(SisterDistance(:,1)),'Color', colormaps_SD(1,:))

        title('Distance to sister cell after division')
        ylabel("µm")
        xlabel("Frames after division")

        hold on
        for i = 2:size(SisterDistance,2)
            plot(1:(length(nonzeros(SisterDistance(:,i)))),nonzeros(SisterDistance(:,i)),'Color', colormaps_SD(i,:))
        end
        hold off
        lgd=legend(num2str(Tracklist(TrackdivisionLogical)),'Location','northeastoutside');
        title(lgd,"Track ID")
    else
        disp('no divisions')
    end
        

    %% plot Neighbour
    
    delta_times = NeighbourPlot(Data,'distall',4,Framelist,Tracklist_sub_name,TrackdivisionLogical,Tracklist);
    NeighbourPlot(Data,'delta',5,Framelist,Tracklist_sub_name,TrackdivisionLogical,Tracklist);
    NeighbourPlot(Data,'del_dvsn',6,Framelist,Tracklist_sub_name,TrackdivisionLogical,Tracklist);
    
    %% calculate a table for elke (cells that have moved a certain distance from start till end point) two tables one for dividing cells and one for non dividing cells.
    %which is also the output
    
    %start by nondividing cells
    Tracklist_nonD = Tracklist(~TrackdivisionLogical);
    moved_non_dividing= zeros(length(Tracklist_nonD),1);
    for j = 1:length(Tracklist_nonD)
       moved_non_dividing(j)= Calculate_moved(Data(Data.TrackID == Tracklist_nonD(j),:),0); 
    end
    [N_nondividing,edges_non]= histcounts(moved_non_dividing,0:5:max(moved_non_dividing));

    %now for dividing 
    %This need to be done better right now it only works for 1 division!!!!
    Tracklist_D = Tracklist(TrackdivisionLogical);
    %Tracklist_D_withsubnotation = unique(Data.TrackID_sub(Data.TrackID == Tracklist_D(1)));
    counter=1;
    moved_dividing = zeros(length(Tracklist_D)*2);
    for l=1:length(Tracklist_D)
        Track_sub = unique(Data.TrackID_sub(Data.TrackID == Tracklist_D(l)));

        if length(Track_sub) > 3
            disp("more than one division is not supported")
        end
        Mother = Data(Data.TrackID_sub == Track_sub(1),:);
        Daughter1 = Data(Data.TrackID_sub == Track_sub(2),:);
        Daughter2 = Data(Data.TrackID_sub == Track_sub(3),:);
        moved_dividing(counter) = Calculate_moved([Mother;Daughter1],0);
        counter = counter+1;
        moved_dividing(counter) = Calculate_moved([Mother;Daughter2],0);
        counter = counter+1;   
    end
    
    if sum(TrackdivisionLogical) ~= 0
        [N_dividing,edges_d]= histcounts(moved_dividing,0:5:max(moved_dividing));
        elkehist = zeros(max([length(edges_non),length(edges_d)]),3);
        elkehist(1:length(edges_d),1)=edges_d;
        elkehist(1:length(N_dividing),3) = N_dividing;
    else
        elkehist = zeros(length(edges_non),3);
    end
    
    %for now asign it out in the worksheet
    
    elkehist(1:length(edges_non),1)=edges_non;
   
    elkehist(1:length(N_nondividing),2) = N_nondividing;
    
    output = elkehist;
end

%% seperate functions:
function output = NeighbourPlot(Data,plottomake,figure_nr,Framelist,Tracklist_sub_name,TrackdivisionLogical,Tracklist)
            
            
    figure(figure_nr);
    clf;
    if ~strcmp(plottomake,'del_dvsn')
        ax_l = subplot(1,3,1);
        ax_m = subplot(1,3,2);
        ax_h = subplot(1,3,3);
    else
        ax_l = subplot(2,3,1);
        ax_m = subplot(2,3,2);
        ax_h = subplot(2,3,3);
        ax_ld = subplot(2,3,4);
        ax_md = subplot(2,3,5);
        ax_hd = subplot(2,3,6);
    end

                
            %fix the axes
            xlim(ax_l, [0, max(Framelist)])
            xlim(ax_m, [0, max(Framelist)])
            xlim(ax_h, [0, max(Framelist)])
            
            if strcmp(plottomake, 'distall')
                ymax = max(Data.FirstNeighbourDistance);
                ymin = 0;
            else
                %need to run through the list to figure out the range
                largestpoint = zeros(length(Tracklist_sub_name),1);
                lowestpoint = zeros(length(Tracklist_sub_name),1);
                for k = 1:length(Tracklist_sub_name)
                    plot_data=[Data.frame(Data.TrackID_sub == Tracklist_sub_name(k)),Data.FirstNeighbourDistance(Data.TrackID_sub == Tracklist_sub_name(k))];
                    largestpoint(k) = max(plot_data(:,2)-plot_data(1,2));
                    lowestpoint(k) = min(plot_data(:,2)-plot_data(1,2));
                end
                ymax = max(largestpoint);
                ymin = min(lowestpoint);
            end
            ylim(ax_l, [ymin, ymax])
            ylim(ax_m, [ymin, ymax])
            ylim(ax_h, [ymin, ymax])
            hold(ax_l,'on')
            hold(ax_m,'on')
            hold(ax_h,'on')
            if strcmp(plottomake,'del_dvsn')
                ylim(ax_ld, [ymin, ymax])
                ylim(ax_md, [ymin, ymax])
                ylim(ax_hd, [ymin, ymax])
                hold(ax_ld,'on')
                hold(ax_md,'on')
                hold(ax_hd,'on')
                lowdivid=0;
                meddivid=0;
                highdivid=0;
            end
            high_threshold = ymax*2/3;
            med_threshold  = ymax*1/3;
            low=0;
            med=0;
            high=0;
            delta_time = zeros(length(Tracklist_sub_name),3);
            for k = 1:length(Tracklist_sub_name)
                plot_data=[Data.frame(Data.TrackID_sub == Tracklist_sub_name(k)),Data.FirstNeighbourDistance(Data.TrackID_sub == Tracklist_sub_name(k))];
            
                if ~strcmp(plottomake,'distall')
                    plot_data(:,2) = plot_data(:,2)-plot_data(1,2);
                end
            
                %determine axis
                if max(plot_data(:,2)) > high_threshold
                    if strcmp(plottomake,'del_dvsn') && TrackdivisionLogical(Tracklist == max(Data.TrackID(Data.TrackID_sub == Tracklist_sub_name(k))))
                        current_axes = ax_hd;
                        highdivid = highdivid +1;
                    else
                        current_axes = ax_h;
                    end
                    high = high +1;
                elseif max(plot_data(:,2)) > med_threshold
                    if strcmp(plottomake,'del_dvsn') && TrackdivisionLogical(Tracklist == max(Data.TrackID(Data.TrackID_sub == Tracklist_sub_name(k))))
                        current_axes = ax_md;
                        meddivid = meddivid + 1;
                    else
                        current_axes = ax_m;
                    end
                    med = med +1;
                else
                    if strcmp(plottomake,'del_dvsn') && TrackdivisionLogical(Tracklist == max(Data.TrackID(Data.TrackID_sub == Tracklist_sub_name(k))))
                        current_axes = ax_ld;
                        lowdivid = lowdivid +1;
                    else
                        current_axes = ax_l;
                    end
                    low = low +1;
                end
                plot(current_axes,plot_data(:,1),plot_data(:,2))
%                 if ~ishold(current_axes)
%                     hold(current_axes,'on')
%                 end
                %delta_time(k,1:2) = [plot_data(end,2)-plot_data(1,2),length(plot_data)];
                nznplotdata= nonzeros(plot_data(:,2));
                if ~isempty(nznplotdata)
                    delta_time(k,1:2) = [nznplotdata(end)-nznplotdata(1),length(plot_data)];
                    delta_time(:,3) = delta_time(:,1)./delta_time(:,2);
                end
            end
            output=delta_time;
            switch plottomake
                case 'distall'
                    name = 'Distance to Neighbour first frame';
                case 'delta'
                    name = 'Distance to Neighbour first frame (normalized)';
                case 'del_dvsn'
                   name = 'Distance to Neighbour first frame (normalized)'; 
            end
            if strcmp(plottomake,'del_dvsn')
                med = med-meddivid;
                high = high-highdivid;
                low = low - lowdivid;
            end
            title(ax_l, {[name];['Close distance (' num2str(low) ' tracks)'];['max distance < ' num2str(med_threshold) ' µm']})
            title(ax_m, {[name];['Medium distance (' num2str(med) ' tracks)']; [num2str(med_threshold) ' µm < max distance < ' num2str(high_threshold) 'µm']})
            title(ax_h, {[name];['High distance (' num2str(high) ' tracks)']; [num2str(high_threshold) ' µm < max distance']})
            xlabel(ax_l,'Frame')
            xlabel(ax_m,'Frame')
            xlabel(ax_h,'Frame')
            ylabel(ax_l,'µm')
            ylabel(ax_m,'µm')
            ylabel(ax_h,'µm')
            if strcmp(plottomake,'del_dvsn')
                title(ax_ld,['Close distance (' num2str(lowdivid) ' tracks)'])
                title(ax_md,['Medium distance (' num2str(meddivid) ' tracks)'])
                title(ax_hd,['High distance (' num2str(highdivid) ' tracks)'])
                ylabel(ax_ld,{'Dividing cells'; 'µm'}, 'fontsize', 12);
                xlabel(ax_ld,'Frame')
                xlabel(ax_md,'Frame')
                xlabel(ax_hd,'Frame')
                ylabel(ax_l,{'non dividing cells'; 'µm'}, 'fontsize', 12)
                ylabel(ax_md,'µm')
                ylabel(ax_hd,'µm') 
            end
end
        
function distance = Calculate_moved(SubSetData,divide_idx)
    first = min(SubSetData.frame);
    last = max(SubSetData.frame);
    if ~divide_idx
        X=SubSetData.pos_x_corr(SubSetData.frame == last);
        X1=SubSetData.pos_x_corr(SubSetData.frame == first);
        Y=SubSetData.pos_y_corr(SubSetData.frame == last);
        Y1=SubSetData.pos_y_corr(SubSetData.frame == first);
        Z=SubSetData.pos_z_corr(SubSetData.frame == last);
        Z1=SubSetData.pos_z_corr(SubSetData.frame == first);
        distance = ((X-X1)^2+(Y-Y1)^2+(Z-Z1)^2)^0.5;
    else
        
    end
end
