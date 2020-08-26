%
% Description: Segment phase-contrast microscopy images manually
% Last Modified: Aug 25, 2020
% Author: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu>
% Notes:
% Microscopy file I/O provided by Bio-Formats toolbox
% Interactive visualization based on 'zoom2cursor' by Brett Shoelson (shoelson@helix.nih.gov)
% Tested on MATLAB R2017b
%

% Requested features (TODO):
% Automatically place vertex points for editing after freehand drawing
% Propagation of segmented boundaries in image stacks

close all; clear all; clc;

global num_segmentations
global stored_masks
global imfreehand_handles

addpath('bfmatlab');

finished_segmentation = false;

while (finished_segmentation == false)

    [imgfile, imgpath] = uigetfile('*.nd2', 'Select a file');
    imgdata = bfopen(strcat(imgpath, imgfile));

    num_series = size(imgdata, 1);
    assert(num_series == 1, 'ERROR: Expected single series acquisition.')

    num_planes = size(imgdata{1, 1}, 1);
    assert(num_planes == 1, 'ERROR: Expected 1 channel in input image.')

    ch1_pixel_data = imgdata{1, 1}{1, 1};
    nd2_img = double(ch1_pixel_data)/double(max(ch1_pixel_data(:)));

    nrows = size(nd2_img, 1);
    ncols = size(nd2_img, 2);
    
    [fpath, fname, fext] = fileparts(strcat(imgpath, imgfile));

    instr_msg = {'[x] Zoom In', '[z] Zoom Out', '[c] Reset Zoom', '[n] New Segmentation'};
    uictrl = figure('Name', 'UI Control Keys', 'NumberTitle', 'off');
    uicontrol(uictrl, 'style', 'text', 'backgroundcolor', [0.9375 0.9375 0.9375],...
            'foregroundcolor', 'k', 'fontsize', 20, 'units', 'normalized',...
            'position', [0.25 0.25 0.5 0.5], 'string', instr_msg);

    num_segmentations = 0;
    stored_masks = {};

    figure;
    imshow(nd2_img);
    get_segmentations

    uiwait(uictrl);

    if num_segmentations == 0

        warning('Did not store any segmentations');

    else
        
        output_folder = strcat(fpath, filesep, fname);
        
        if exist(output_folder, 'dir')
            rmdir(output_folder, 's')
        end
        
        mkdir(output_folder)
        
        fig2 = figure(2);
        imshow(nd2_img);
        hold on

        output_fileID = fopen(strcat(output_folder, filesep, 'measurements.csv'), 'w');
        fprintf(output_fileID, 'Cell_ID, Centroid_X, Centroid_Y, Area_Px\n');
        output_img = zeros(nrows, ncols, 3);
        
        cell_id = 1;
        bin_mask = stored_masks{1};
        measurements = regionprops(bin_mask);
        assert(size(measurements, 1) == 1)
        
        save(strcat(output_folder, filesep, 'CellID_', int2str(cell_id), '.mat'), 'bin_mask');
        annotations = {int2str(cell_id)};
        centroid_coords = measurements.Centroid;
        fprintf(output_fileID, '%d, %f, %f, %d\n', cell_id, centroid_coords(1), centroid_coords(2), measurements.Area);
        
        cell_bdy_objs = bwboundaries(bin_mask);
        assert(length(cell_bdy_objs) == 1);
        cell_bdy = cell_bdy_objs{1};
        plot(cell_bdy(:,2), cell_bdy(:,1), 'w', 'LineWidth', 2)

        if num_segmentations >= 2
            for i = 2:num_segmentations

                cell_id = cell_id + 1;
                bin_mask = stored_masks{i};
                measurements = regionprops(bin_mask);
                assert(size(measurements, 1) == 1)
                
                save(strcat(output_folder, filesep, 'CellID_', int2str(cell_id), '.mat'), 'bin_mask');
                annotations(end+1) = {int2str(cell_id)};
                cell_centroid = measurements.Centroid;
                centroid_coords = [centroid_coords; cell_centroid];
                fprintf(output_fileID, '%d, %f, %f, %d\n', cell_id, centroid_coords(1), centroid_coords(2), measurements.Area);
                
                cell_bdy_objs = bwboundaries(bin_mask);
                assert(length(cell_bdy_objs) == 1);
                cell_bdy = cell_bdy_objs{1};
                plot(cell_bdy(:,2), cell_bdy(:,1), 'w', 'LineWidth', 2)

            end
        end

        scatter(centroid_coords(:,1), centroid_coords(:,2), 'r', 'filled')
        txt_disp_offset = 10;
        text(centroid_coords(:,1)+txt_disp_offset, centroid_coords(:,2)-txt_disp_offset, annotations, 'Color', 'magenta', 'FontSize', 16)
        hold off
        saveas(gcf, strcat(output_folder, filesep, 'annotated_mask.png'));

    end
    
    usr_response = questdlg('Would you like to open another image?', 'Segmentation Complete', 'Yes', 'No', 'No');
    
    switch usr_response
        case 'No'
            close(fig2);
            finished_segmentation = true;
        case 'Yes'
            close(fig2);
    end
    
    fclose(output_fileID);
    
end

function get_segmentations

    currfig = findobj('type', 'figure');
    
    if isempty(currfig)
        beep;
        h = warndlg('There are no figures open!');
        uiwait(h);
        return
    end
    
    currfig = currfig(1);
    figure(currfig);
    
    zoomparams.currax = get(currfig, 'currentaxes');
    if isempty(zoomparams.currax)
        error('The current figure contains no axes!');
    end
    
    % Precedence: Images, surfaces, lines, patches
    zoomparams.currobj = findobj(currfig, 'type', 'image');
    if isempty(zoomparams.currobj)
        zoomparams.currobj = findobj(currfig, 'type', 'surface');
    end
    if isempty(zoomparams.currobj)
        zoomparams.currobj = findobj(currfig, 'type', 'line');
    end
    if isempty(zoomparams.currobj)
        zoomparams.currobj = findobj(currfig, 'type', 'patch');
    end
    if isempty(zoomparams.currobj)
        error('The current axis doesn''t appear to contain any valid images, surfaces, lines, or patches.');
    end
    
    % Default zoom value is 50%
    zoomparams.pct = 0.5;
    
    zoomparams.currobj = zoomparams.currobj(1);
    zoomparams.objtype = get(zoomparams.currobj, 'type');
    zoomparams.bdfcnold = get(currfig, 'keypressfcn');
    zoomparams.baold = get(zoomparams.currobj, 'busyaction');
    zoomparams.oldpointer = get(currfig, 'pointer');
    
    set(currfig, 'Pointer', 'crosshair');
    axes(zoomparams.currax);
    currax = get(zoomparams.currax);
    zoomparams.refax = copyobj(zoomparams.currax, gcf);
    
    % Ignore "Unrecognized OpenGL" message
    warning off;
    
    set(currfig, 'keypressfcn', 'feval(getappdata(gcf,''keypressfcn''));', 'busyaction', 'queue');
    set(findobj(zoomparams.refax, 'type', 'children'), 'handlevisibility', 'on');
    set(zoomparams.refax, 'visible', 'off');
    axes(zoomparams.refax);
    cla;
    
    zoomparams.oldaxunits = get(zoomparams.currax, 'units');
    zoomparams.ydir = get(zoomparams.currax, 'ydir');
    zoomparams.oldxlim = get(zoomparams.currax, 'xlim');
    zoomparams.oldylim = get(zoomparams.currax, 'ylim');
    zoomparams.oldzlim = get(zoomparams.currax, 'zlim');
    zoomparams.dbold = get(currfig, 'doublebuffer');
    zoomparams.xrange = diff(zoomparams.oldxlim);
    zoomparams.yrange = diff(zoomparams.oldylim);
    zoomparams.zrange = diff(zoomparams.oldzlim);
    zoomparams.xdist = zoomparams.pct*zoomparams.xrange;
    zoomparams.ydist = zoomparams.pct*zoomparams.yrange;
    zoomparams.zdist = zoomparams.pct*zoomparams.zrange;
    zoomparams.oldwbmf = get(currfig, 'windowbuttonmotionfcn');
    
    setappdata(currfig, 'zoomfcnhandle', @zoomfcn);
    setappdata(currfig, 'keypressfcn', @keyfcn);
    
    set(currfig, 'doublebuffer', 'on', 'windowbuttonmotionfcn', 'feval(getappdata(gcf,''zoomfcnhandle''));');
    
    endbutton = uicontrol('style', 'pushbutton', 'string', 'Close','units', 'normalized', 'position', [0 0 0.06 0.045],...
        'foregroundcolor', 'k', 'backgroundcolor', 'w', 'fontsize', 13,...
        'fontweight', 'b', 'callback',...
        ['zoomparams = getappdata(gcf,''zoomparams'');set(gcf,''windowbuttonmotionfcn'',zoomparams.oldwbmf);',...
            'set(zoomparams.currax,''units'',zoomparams.oldaxunits,''xlim'',zoomparams.oldxlim,''ylim'',zoomparams.oldylim);',...
            'set(get(zoomparams.currax,''children''),''buttondownfcn'',zoomparams.bdfcnold,''busyaction'',zoomparams.baold);',...
            'set(gcf,''pointer'',zoomparams.oldpointer,''doublebuffer'',zoomparams.dbold);,',...
            'delete(zoomparams.dispbox1);delete(zoomparams.dispbox2);delete(gcbo);delete(uictrl);']);
        
    zoomparams.dispbox1 = uicontrol('style', 'frame', 'backgroundcolor', 'w', 'units', 'normalized', 'position', [0.0475 -1.5 0.35 0.065]);
    
    if ~strcmp(zoomparams.objtype, 'surface')
        msgstr = sprintf('Cursor X = %3.0f  Cursor Y = %3.0f', 0, 0);
    else
        msgstr = sprintf('Cursor X = %3.0f  Cursor Y = %3.0f Z = %3.0f', 0, 0, 0);
    end
    
    zoomparams.dispbox2 = uicontrol('style', 'text', 'backgroundcolor', [0.9375 0.9375 0.9375],...
        'foregroundcolor', 'k', 'units', 'normalized',...
        'position', [0.065 0.0005 0.32 0.04], 'fontsize', 14, 'string', msgstr,...
        'horizontalalignment', 'l');
    
    setappdata(gcf, 'zoomparams', zoomparams);
    
end
    
function zoomfcn

    zoomparams = getappdata(gcf, 'zoomparams');
    posn = get(zoomparams.refax, 'currentpoint');
    posn = posn(1,:);
    x = posn(1,1);
    y = posn(1,2);
    z = posn(1,3);
    
    switch zoomparams.objtype
        
        case 'image'
            
            % x and y are already in expressed in proper pixel coordinates
            x1 = min(max(1,x-0.5*zoomparams.xdist), zoomparams.xrange-zoomparams.xdist) + 0.5;
            y1 = min(max(1,y-0.5*zoomparams.ydist), zoomparams.yrange-zoomparams.ydist) + 0.5;
            z1 = min(max(1,z-0.5*zoomparams.zdist), zoomparams.zrange-zoomparams.zdist) + 0.5;
            x2 = x1 + zoomparams.xdist;
            y2 = y1 + zoomparams.ydist;
            z2 = z1 + zoomparams.zdist;
            
        case {'line', 'surface', 'patch'}
            
            % x, y and z are in normalized units; must be converted
            x = zoomparams.oldxlim(1) + x*zoomparams.xrange;
            y = zoomparams.oldylim(1) + y*zoomparams.yrange;
            z = zoomparams.oldzlim(1) + z*zoomparams.zrange;
       
            % now change the limits of currax, ensuring that the original limits are not exceeded
            x1 = max(x-zoomparams.xdist/2,zoomparams.oldxlim(1));
            y1 = max(y-zoomparams.ydist/2,zoomparams.oldylim(1));
            z1 = max(z-zoomparams.zdist/2,zoomparams.oldzlim(1));
            x2 = x1+zoomparams.xdist;
            y2 = y1+zoomparams.ydist;
            z2 = z1+zoomparams.zdist;
            
            % if new limits are out of range, adjust them:
            if x2 > zoomparams.oldxlim(2)
                x2 = zoomparams.oldxlim(2);
                x1 = x2 - zoomparams.xdist;
            end
            if y2 > zoomparams.oldylim(2)
                y2 = zoomparams.oldylim(2);
                y1 = y2 - zoomparams.ydist;
            end
            if z2 > zoomparams.oldzlim(2)
                z2 = zoomparams.oldzlim(2);
                z1 = z2 - zoomparams.zdist;
            end
            
            % now get the x,y positions in currax for display purposes
            posn = get(zoomparams.currax, 'currentpoint');
            posn = posn(1,:);
            x = posn(1,1);
            y = posn(1,2);
            z = posn(1,3);
    end
    
    if x >= zoomparams.oldxlim(1) && x <= zoomparams.oldxlim(2) && ...
            y >= zoomparams.oldylim(1) && y <= zoomparams.oldylim(2) && ...
        z >= zoomparams.oldzlim(1) && z <= zoomparams.oldzlim(2)
        if strcmp(zoomparams.objtype, 'surface')
            set(zoomparams.dispbox2, 'string', sprintf('Cursor X = %3.2f  Cursor Y = %3.2f Z = %3.2f', x, y, z));
            set(zoomparams.currax, 'xlim', [x1 x2], 'ylim', [y1 y2], 'zlim', [z1 z2]);
        else
            set(zoomparams.dispbox2, 'string', sprintf('Cursor X = %3.2f  Cursor Y = %3.2f', x, y));
            set(zoomparams.currax, 'xlim', [x1 x2], 'ylim', [y1 y2]);
        end
    else
        if strcmp(zoomparams.objtype, 'surface')
            set(zoomparams.dispbox2, 'string', sprintf('Cursor X = %3.f  Cursor Y = %3.0f Z = %3.0f', 0, 0, 0));
            set(zoomparams.currax, 'xlim', zoomparams.oldxlim, 'ylim', zoomparams.oldylim, 'zlim', zoomparams.oldzlim);
        else
            set(zoomparams.dispbox2, 'string', sprintf('Cursor X = %3.0f  Cursor Y = %3.0f', 0, 0));
            set(zoomparams.currax, 'xlim', zoomparams.oldxlim, 'ylim', zoomparams.oldylim);
        end
    end
    
    return
    
end

function keyfcn

    global num_segmentations
    global stored_masks
    global imfreehand_handles

    zoomparams = getappdata(gcf, 'zoomparams');
    key_pressed = get(gcf, 'CurrentKey');
    
    switch key_pressed
        
        case 'z'
            
            zoomparams.pct = min(1, zoomparams.pct*1.1);
            
        case 'x'
            
            zoomparams.pct = max(0.01, zoomparams.pct*0.9);
            
        case 'c'
            
            zoomparams.pct = 1;
            
        case 'd'
            
            assert(length(imfreehand_handles) == length(stored_masks));
            assert(length(stored_masks) == num_segmentations)
            
            if ~isempty(imfreehand_handles)
                
                bdy =  imfreehand_handles(end);
                delete(bdy);
                imfreehand_handles(end) = [];
                stored_masks(end) = [];
                num_segmentations = num_segmentations - 1;
                
            end
        
        case 'n'
            
            setappdata(gcf, 'zoomfcnhandle', @() drawnow('update'));
            bdy = imfreehand(zoomparams.currax);
            bdy.Deletable = false;
            bin_mask = createMask(bdy, zoomparams.currobj);
            imfreehand_handles = [imfreehand_handles bdy];
            num_segmentations = num_segmentations + 1;
            stored_masks{num_segmentations} = bin_mask;
            setappdata(gcf, 'zoomfcnhandle', @zoomfcn);

    end
    
    zoomparams.xdist = zoomparams.pct*zoomparams.xrange;
    zoomparams.ydist = zoomparams.pct*zoomparams.yrange;
    zoomparams.zdist = zoomparams.pct*zoomparams.zrange;
    if ~strcmp(get(gcf, 'selectiontype'), 'extend')
        setappdata(gcf, 'zoomparams', zoomparams);
        feval(getappdata(gcf, 'zoomfcnhandle'));
    end
    
    return
    
end