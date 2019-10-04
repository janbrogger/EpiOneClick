# EpiOneClick
A method to annotate epileptiform elements with one click.
Developed by Eivind Aanestad
Not yet completed for easy use by others, but get in touch if you're interested

How to set up epileptiform annotation

Prerequisities:
1. MATLAB
2. An SQL Server

Procedure
0. Checkout repos
  cd "C:\Midlertidig Lagring\"
  git clone https://github.com/janbrogger/EpiOneClick epileptiform
  git clone https://github.com/janbrogger/ScorePipeline ScorePipeline
  git clone https://github.com/sccn/eeglab.git eeglab
  git clone https://github.com/janbrogger/fieldtrip fieldtrip
  
2. Setup the ScorePipeline according to the instructions in its setup.docx

a) take a SCORE database, mount it

b) run the "create tables" scripts from ScorePipeline

e) Run the script ScorePipeline\update-fieldtrip-in-eeglab to update fieldtrip in EEGLAB

d) Start Matlab and enter the command: cd [path-to-ScorePipeline\matlab]

e) Then the command

ScorePipeline

f) You probably have to run the command ScoreAssistJDBC and restart matlab

g) Then ScorePipeline again

3. The epileptiform annotation should have already been set up.

In SQL Server management studio, in the table SearchResult, field MouseDownFunction

should contain the contents of the file ScoreMouseDownForPaper2.m

