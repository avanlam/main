function codes = colour_groups(shade, n)
%COLOUR_GROUPS : Returns gradients of n colours for the given shade, as an
% array of nx3

% Pantone colours for Spring/Summer 2021
BurntCoral = [233,137,126]./255;
Raspberry = [210,56,108]./255;
Rust = [181, 90,48]./255;

Illuminating = [245,223,77]./255;
Marigold = [253,172,83]./255;

GreenAsh = [160, 218, 169]./255;
Mint = [0,161,112]./255;

Cerulean = [155,183,212]./255;
FrenchBlue = [0,114,181]./255;

Lavender = [226 172 215]./255;
Orchid = [146,106,166]./255;

if shade == "blue"
    codes = [linspace(150, 0, n)', linspace(180, 120, n)', linspace(220, 180, n)'];
elseif shade == "green"
    codes = [linspace(160, 0, n)', linspace(220, 160, n)', linspace(170, 110, n)'];
elseif shade == "red"
    codes = [linspace(230, 210, n)', linspace(135, 60, n)', linspace(130, 110, n)'];
elseif shade == "orange"
    codes = [linspace(240, 255, n)', linspace(200, 130, n)', linspace(100, 80, n)'];
elseif shade == "purple"
    codes = [linspace(220, 145, n)', linspace(180, 110, n)', linspace(220, 170, n)'];
end
        
codes = codes./255;
end

