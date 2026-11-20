function app = Start_HX711_BPA()
    %START_HX711_BPA Launch the BPA HX711 force/pressure test app.
    %
    % Works from any current folder and any clone location on any machine:
    % it adds THIS file's folder to the MATLAB path (which also makes the
    % Arduino add-on package +arduinoioaddons/+basicHX711 resolvable for
    % arduino(...,'libraries',{'basicHX711/basic_HX711'})), then starts the
    % app. Because the app class is uniquely named HX711_BPA, it can never
    % collide with any installed File-Exchange "HX711" app or add-on copy.
    %
    % Prerequisites on each machine:
    %   - MATLAB R2025a or newer
    %   - "MATLAB Support Package for Arduino Hardware" installed
    %     (Add-Ons > Get Hardware Support Packages)
    %   - the HX711 Arduino library is uploaded automatically from this
    %     folder on first Connect (no manual Arduino IDE step needed).

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    app = HX711_BPA(here);
    if nargout == 0
        clear app
    end
end
