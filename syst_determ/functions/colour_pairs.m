function codes = colour_pairs(shade)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Rust = [181, 90,48]./255;
BurntCoral = [233,137,126]./255;

Raspberry = [210,56,108]./255;

Marigold = [253,172,83]./255;
Illuminating = [245,223,77]./255;

GreenAsh = [160, 218, 169]./255;
Mint = [0,161,112]./255;

Cerulean = [155,183,212]./255;
FrenchBlue = [0,114,181]./255;

Lavender = [226 172 215]./255;
Orchid = [146,106,166]./255;

if shade == "blue"
    codes{1} = Cerulean;
    codes{2} = FrenchBlue;
elseif shade == "green"
    codes{1} = GreenAsh;
    codes{2} = Mint;
elseif shade == "red"
    codes{1} = [251,154,153]./255;
    codes{2} = [227,26,28]./255;
elseif shade == "orange"
    codes{1} = Illuminating;
    codes{2} = Marigold;
elseif shade == "purple"
    codes{1} = Lavender;
    codes{2} = Orchid;
elseif shade == "brown"
    codes{1} = BurntCoral;
    codes{2} = Rust;
elseif shade == "spring"
    codes{1} = GreenAsh;
    codes{2} = Mint;
    codes{3} = Marigold;
elseif shade == "fresh"
    codes{1} = GreenAsh;
    codes{2} = Illuminating;
    codes{3} = Marigold;
elseif shade == "duo"
    codes{1} = Marigold;
    codes{2} = Mint;
elseif shade == "all"
    
    codes{1} = BurntCoral;
    codes{2} = Marigold;
    codes{3} = Illuminating;
    codes{4} = GreenAsh;
    codes{5} = Mint;
    codes{6} = Cerulean;
    codes{7} = FrenchBlue;
    codes{8} = Lavender;
    codes{9} = Orchid;
    
elseif shade == "three"
    codes{1} = GreenAsh;
    codes{2} = Mint;
    codes{1} = Illuminating;
end
        
end

