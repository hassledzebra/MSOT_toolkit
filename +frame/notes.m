clear all; close all; clc;
% It would be relevant for there to be a Coordinate object, with a 'natural'
% interpolation scheme attached. This could also be related to the 'natural'
% coordinate system for that object, and could be as complicated as e.g.
% 6-sphere coordinates or Plucker coordinates. This would likely require a
% registry of names and interconversions. 


%% For transducers themselves (index-based)

Nxdcr = 256;
xdcr_index = 1:Nxdcr;
xdcr_R = 0.04055;
    xdcr_R = xdcr_R(:) .* ones(Nxdcr,1); % Singleton expansion of coordinate.
xdcr_th = linspace(3*pi/4,9*pi/4,256);

xdcr_interp = 'cylindrical';
Coord_xdcr_Dim_Names = {'CYLINDRICAL_R', 'CYLINDRICAL_THETA'}; % Note standardization
Coord_xdcr_Dim_Units = {'METERS','RADIANS'}; % Note standardization
Coord_xdcr_Name = {'TRANSDUCER'};
Coord_xdcr_relationship = {'CYLINDRICAL_R(+)CYLINDRICAL_THETA'}; %(+) denotes a direct sum, where the coordinates are tacked together pairwise.
                                                                % This seems to
                                                                % not strictly
                                                                % be a direct
                                                                % sum, but I
                                                                % don't know
                                                                % what else to
                                                                % call it.
                                                                % The main
                                                                % requirement
                                                                % here is that
                                                                % both
                                                                % coordinates
                                                                % have the same
                                                                % dimension.

N_xdcr_coords = size(Coord_xdcr_Dim_Names,2); 

xdcr_coords = {xdcr_R(:),xdcr_th(:)}; % 'Coordinate-first' definition - this is a 
                                % 1 x 2 cell array, made by just gluing together
                                % the lists of coordinates.

xdcr_index_coords = mat2cell(cell2mat(xdcr_coords) , ones(Nxdcr,1) , N_xdcr_coords); % 'Index-first' 
                                % definition - this is a 256 x 1 cell array,
                                % made by gluing together the coordinates for
                                % each transducer individually. For transducers,
                                % this is probably the more 'natural'
                                % definition.
                                
%% For the sampling times

Nt = 2030;

t_index = 1:Nt;
dt = 1/40E6;
t = t_index.*dt;

t_interp = 'linear';
Coord_Dim_Names = {'TIME'}; % Note standardization
Coord_Dim_Units = {'SECONDS'}; % Note standardization
Coord_Name = {'SAMPLING_TIME'}; % Free to determine.
N_t_coords = size(Coord_Dim_Names,2); 

t_coords = {t(:)};
t_index_coords = mat2cell(cell2mat(t_coords), ones(Nt,1), N_t_coords);


%% For the data frame itself.


Coord_Dim_Names = {'SAMPLING_TIME','TRANSDUCER'};
Coord_relationship = {'SAMPLING_TIME(x)TRANSDUCER'}; % (x) denotes the tensor product, where we're taking all pairs of coords.

Coord_size = [Nt,Nxdcr];
N_sample_coords = prod(Coord_size); % prod because tensor product makes a space with product size.
sample_index = 1:N_sample_coords;

sample_coords = {t_coords,xdcr_coords};
sample_coords_2 = {cell2mat(t_coords),cell2mat(xdcr_coords)};

[T_grid,Xdcr_grid] = ndgrid(1:Nt,1:Nxdcr);

coords_for_each_sample = [sample_coords_2{1}(T_grid(:),:) sample_coords_2{2}(Xdcr_grid(:),:)];








                               
                                

