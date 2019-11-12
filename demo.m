settings = OneClickSettings();
addpath(settings.EeglabPath)
OpenEegFileInEeglab('janbrogger.e');
GotoEvent(40);
