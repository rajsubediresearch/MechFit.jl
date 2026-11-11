
function adjustedFontSize=GetAdjustedFontSize()

% === Retrieve screen DPI and size ===
screenDPI = get(0, 'ScreenPixelsPerInch');
screenSize = get(0, 'ScreenSize');  % [left, bottom, width, height]

% === Define base values ===
baseDPI = 96;              % Standard DPI for which baseFontSize is designed
baseFontSize = 24;         % Ideal font size for base DPI
maxFontSize = 24;          % Upper limit to avoid oversized fonts
minFontSize = 10;          % Lower limit for legibility

% === Scale font size according to screen DPI ===
scalingFactor = screenDPI / baseDPI;
adjustedFontSize = round(baseFontSize * scalingFactor);

% === Clamp the font size to a reasonable range ===
adjustedFontSize = min(max(adjustedFontSize, minFontSize), maxFontSize);