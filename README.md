# EpiOneClick
A method to annotate epileptiform elements with one click.
Developed by Eivind Aanestad
Not yet completed for easy use by others, but get in touch if you're interested

Take a look at the demo.m script.

# Instructions for use 
Look up the "demo" script.

1. Change the Matlab path to the directory of this tool
2. Edit folders in "EpiOneClickSettings.m".
3. Call the EpiOneClick method, for example: EpiOneClick('example-eeg\janbrogger.e');

# Bullet points for successful annotation:
-Choose optimal channel.

-Choose optimal wave.
Click on the peak of the waveform you consider to be the most epileptiform with the criteria in mind.

-Check result.
A window will pop up showing the result of auto-annotation. If the spike peak hits the mark, the annotation was successful. If not, close pop-up window and try again. Click more precisely close to the peak.

-Autoannotation pop-up-window legend:
First blue mark: Spike start. 
Red mark: Peak. 
Second blue mark: Spike end and slow-wave start. 
Third blue mark: Slow-wave end.
Additionally: The green dots show other local minima considered as waveform limits, a green line showing a smoothed helping signal to find slow-wave end, and a red line showing the Gauss-fit.
