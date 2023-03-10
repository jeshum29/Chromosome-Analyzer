function writeTIFFmeta(dd, MIN_INTENSITY, MAX_INTENSITY, Z_RES)

% dd is the name of a tiff file

if ~exist('Z_RES', 'var')
    MIN_INTENSITY = 1;
    MAX_INTENSITY = 3333;
    Z_RES = 3;
end


minStr = num2str(MIN_INTENSITY);
maxStr = num2str(MAX_INTENSITY);
zResStr = num2str(Z_RES);

t = Tiff(dd, 'r+');

for ii = 1:9999
    
    try
        t.setDirectory(ii);
    catch
        numDir = ii - 1;
        break;
    end
    
end

numDirStr = num2str(numDir);
slashn = '\n';

str = [...
    'ImageJ=1.45b' slashn ...
    'images=' numDirStr slashn...
    'slices=' numDirStr slashn...
    'unit=pixel' slashn...
    'spacing=' zResStr slashn...
    'loop=false' slashn...
    'min=' minStr slashn...
    'max=' maxStr slashn...
    ];

tagStruct.ImageDescription = sprintf(str);


for ii = 1:numDir
    
    t.setDirectory(ii);
    t.setTag(tagStruct);
    t.rewriteDirectory();

end

t.close();

