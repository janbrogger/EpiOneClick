settings = OneClickSettings();
addpath(settings.EeglabPath)
%OpenEegFileInEeglab('janbrogger.e');
OpenEegFileInEeglab('\\hbemta-nevrofil01.knf.local\Rutine\workarea\74399\Patient65_EEG224-OLD_t1.e');
GotoEvent(40);
