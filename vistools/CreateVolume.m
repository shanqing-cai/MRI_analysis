function CreateVolume(brainVolFN, coordFN, outputVolFN, varargin)
%% Config
DEFAULT_SHAPE = 's';
DEFAULT_SIZE = 11;
DEFAULT_CROSSHAIR_THICKNESS = 3;

%% Check dependencies
rpath = which('MRIread');
if isempty(rpath)
    error('Cannot find path to required program: MRIread');
end

rpath = which('MRIwrite');
if isempty(rpath)
    error('Cannot find path to required program: MRIwrite');
end

%% Check input arguments
if nargin <  3
    fprintf(1, 'USAGE: CreateVolume(brainVolFN, coordFN, outputVolFN);\n');
    fprintf(1, '       CreateVolume(brainVolFN, coordFN, outputVolFN, shape);\n');
    fprintf(1, '       CreateVolume(brainVolFN, coordFN, outputVolFN, shape, size);\n');
    fprintf(1, '\nNOTE: In the 2nd and 3rd usage, the shape and/or size will be overridden\n');
    fprintf(1, '        by entries in the coodinates file.\n');
    fprintf(1, '\nshape:\n');
    fprintf(1, '    s - sphere (size --> side length, in mm)\n');
    fprintf(1, '    c - crosshair (size --> length of each axis, in mm\n');
    fprintf(1, '    u - cube (size --> diameter, in mm)\n');
    return
end

if nargin > 3
    DEFAULT_SHAPE = varargin{1};
    
    if nargin > 4
        DEFAULT_SIZE = varargin{2};
    end
end

%% Read coordinates file
if isequal(coordFN(end - 3 : end), '.txt')
    crdInfo = readCoordsFile(coordFN, DEFAULT_SHAPE, DEFAULT_SIZE);
elseif isequal(coordFN(end - 3 : end), '.xls')
    crdInfo = readCoordsFile_xls(coordFN, DEFAULT_SHAPE, DEFAULT_SIZE);
else
    error('Unrecognized file extension name: %s', crdInfo(end - 3 : end));
end


assert(size(crdInfo.crd, 1) == length(crdInfo.val));
assert(size(crdInfo.crd, 1) == length(crdInfo.shape));
assert(size(crdInfo.crd, 1) == length(crdInfo.size));

if size(crdInfo.crd, 1) == 0
    error('Coordinate file %s contains no points', coordFN);
end

%% Read brain volume
imb = MRIread(brainVolFN);

if length(unique(imb.volres)) ~= 1
    error('The input brain volume has non-isotropic voxels.');
end

vs = imb.volres(1);

%% Filling
imc = imb;
imc.vol(:) = 0;
for i1 = 1 : length(crdInfo.val)
    t_mni = crdInfo.crd(i1, :);
    t_idx = round(mni2raster(t_mni));
    
    siz = crdInfo.size(i1) / vs; % Unit: vox
    
    ix = t_idx(1) - floor(siz / 2) : t_idx(1) - floor(siz / 2) + siz - 1;
    iy = t_idx(2) - floor(siz / 2) : t_idx(2) - floor(siz / 2) + siz - 1;
    iz = t_idx(3) - floor(siz / 2) : t_idx(3) - floor(siz / 2) + siz - 1;
    
    if isequal(crdInfo.shape{i1}, 'u') % -- Cube -- %            
        imc.vol(ix, iy, iz) = crdInfo.val(i1);
        
    elseif isequal(crdInfo.shape{i1}, 'c') % -- 3D crosshair -- %        
        wh = (DEFAULT_CROSSHAIR_THICKNESS - 1) / 2;
        
        imc.vol(ix, t_idx(2) + [-wh : wh], t_idx(3) + [-wh : wh]) = crdInfo.val(i1);
        imc.vol(t_idx(1) + [-wh : wh], iy, t_idx(3) + [-wh : wh]) = crdInfo.val(i1);
        imc.vol(t_idx(1) + [-wh : wh], t_idx(2) + [-wh : wh], iz) = crdInfo.val(i1);
        
    elseif isequal(crdInfo.shape{i1}, 's') % -- Shere -- %
        for kx = ix(1) : ix(end)
            for ky = iy(1) : iy(end)
                for kz = iz(1) : iz(end)
                    if (kx - t_idx(1)) ^ 2 ...
                       + (ky - t_idx(2)) ^ 2 ...
                       + (kz - t_idx(3)) ^ 2 <= (siz / 2) ^ 2
                        imc.vol(kx, ky, kz) = crdInfo.val(i1);
                    end
                end
            end
        end
        
    end
end

%% Write to output volume file
if exist(outputVolFN) == 2
    delete(outputVolFN);
end

MRIwrite(imc, outputVolFN);

return

%% Sub-routines
function coord = mni2raster(mniCoord)
assert(length(mniCoord) == 3);

coord = nan(1, 3);
coord(2) = 127.0 - mniCoord(1);
coord(1) = 147.0 - mniCoord(3);
coord(3) = 145.0 + mniCoord(2);
return

%%
function crdInfo = readCoordsFile(coordFN, defShape, defSize)
crdInfo = struct();
crdInfo.crd = nan(0, 3);
crdInfo.val = [];
crdInfo.shape = {};
crdInfo.size = [];

ctxt = textread(coordFN, '%s', 'delimiter', '\n');

for i1 = 1 : numel(ctxt)
    cline = deblank(ctxt{i1});
    cline0 = cline;
    if isempty(cline)
        continue;
    end
    
    if isequal(cline(1), '#')
        continue;        
    end
    
    if ~isempty(strfind(cline, '#'))
        citems = splitstring(cline, '#');
        cline = citems{1};
    end
    
%     ncom = length(strfind(cline, ','));
%     if ncom < 3 || ncom > 5
%         error('Unrecognized format in coordinates file line "%s"', cline0);
%     end
    
    cline = strrep(cline, ',', ' ');
    c_items = splitstring(cline, ' ');
    
    if numel(c_items) < 4 || numel(c_items) > 6
        error('Unrecognized format in coordinates file line "%s"', cline0);
    end
        
    t_coord = [str2double(deblank(c_items{1})), ...
               str2double(deblank(c_items{2})), ...
               str2double(deblank(c_items{3}))];
	t_val = str2double(deblank(c_items{4}));
    
    if t_val == 0
        fprintf(1, 'WARNING: value == 0 for coordinate [%.1f, %.1f, %.1f]', ...
                t_coord(1), t_coord(2), t_coord(3));
    end
           
	if numel(c_items) > 4
        t_shape = strrep(strrep(deblank(c_items{5}), ' ', ''),  '\t', '');
        
        if numel(c_items) > 5
            t_size = str2double(deblank(c_items{6}));
        else
            t_size = defSize;
        end
    else
        t_shape = defShape;
        t_size = defSize;
    end
    
    crdInfo.crd = [crdInfo.crd; t_coord];
    crdInfo.val(end + 1) = t_val;
    crdInfo.shape{end + 1} = t_shape;
    crdInfo.size(end + 1) = t_size;
end
return

%% 
function crdInfo = readCoordsFile_xls(coordFN, defShape, defSize)
crdInfo = struct();
crdInfo.crd = nan(0, 3);
crdInfo.val = [];
crdInfo.shape = {};
crdInfo.size = [];

[N, T] = xlsread(coordFN);

for i1 = 1 : size(N, 1)
    nline = N(i1, :);
    if i1 + 1 <= size(T, 1);
        tline = T(i1 + 1, :);
    else
        tline = {};
    end
    
    if length(nline) < 4
        error('Input xls file format error');
    end
        
    t_coord = [nline(1), nline(2), nline(3)];
	t_val = nline(4);
    
    if t_val == 0
        fprintf(1, 'WARNING: value == 0 for coordinate [%.1f, %.1f, %.1f]', ...
                t_coord(1), t_coord(2), t_coord(3));
    end
           
	if numel(tline) > 4
        t_shape = tline{5};        
        
        if numel(nline) > 5 && ~isnan(nline(6))
            t_size = nline(6);
        else
            t_size = defSize;
        end
    else
        t_shape = defShape;
        t_size = defSize;
    end
    
    crdInfo.crd = [crdInfo.crd; t_coord];
    crdInfo.val(end + 1) = t_val;
    crdInfo.shape{end + 1} = t_shape;
    crdInfo.size(end + 1) = t_size;
end
return