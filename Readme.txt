
(1) Citation to the manuscript: This document describes the data files and code that accompany ``Long-term cycles in food producing archaeological regions" by Jacob Freeman, John M. Anderies, and Jacopo A. Baggio. This README document describes the csv files and how to get started with the analysis. 


(2) Description of the files included in the repository and their relationship to the figures and tables in the manuscript.

Directories:


(i) Modelcode directory: These files are needed to analyze the model numerically and replicate figures 1-3 in the main manuscript. Claude was used to assist with the .jl files to construct figs 2 and 3. CLAUDE.md provides meta data for the construction of the figures. 


(ii)RadiocarbonDataCode directory:

These files contains all radiocarbon ages associated with archaeological remains collected from the archaeological regions. These files are essential to reproduce each file called RegionPerCap.csv and, ultimately, Fig. 4 in the main text. These data are provided for researchers interested in reproducing the KDEs used in the study and engaging in their own analysis. To construct KDEs of radiocarbon ages in each region, researchers should follow the code in CaseStudyRadiocarbon26.R. The raw radiocarbon data are found in RawNA1.scv and RawP3kc14.scsv.

Typical Columns in the radiocarbon files include:

Latitude-decimal degrees
Longitude-decimal degrees
SiteName--name of archaeological site
SiteID" Unique id number for an archaeological site.	
Trinomial--unique trinomial for each site	
Assay No.--unique id for lab and each radiocarbon sample	
Provenience--provenience of dated material if known
Feature Type--archaeological feature from which material was recovered
Material--material dated	
Taxa dated-species or genera of data material	
Human--yes or no. yes indicates human remains dated	
Age-- radiocarbon age	
Error--standard deviation of raw radiocarbon age	
Corrected/Raw	
Corrected/Normalized 14C age	 
Corrected/Normalized Age	
Delta 13C	
Delta 13C source	
Raw/Measured Age	 
Raw/Measured Age	
AMS or Radiometric	
Comments	
Reference


 File path: RadiocarbonDataCode/data/KDEs/`region'.csv 

These files contains the time-series of all simulated KDEs, the mean KDE, and per capita growth rate of the mean KDE. The files are essential to reproduce the economic output dynamics of Figure 4 in the main manuscript and the figures in the Supporting Material. These files can be used by those people who do not want to reproduce the full analysis of archaeological radiocarbon based on the raw data.

Each regional file contains the same columns. For example, NAKDE50Bin.csv contains the following columns for the KDEs in the continental US:

calBP--30 year time-steps in cal BP.
PeriodID--Cultural historical time period	
v1-v500: The summed values of each KDE in 30 year intervals
MKDE--The mean value of the 500 KDEs at each 30 year time-step
PerCap--The per capita growth rate of the mean KDE


(iii) ModelCalibration directory: These files contain python code to calibrate the model to the mean KDE of each archaeological region. As noted in the main manuscript, we provide code to fit a ``null model" in which the parameters of the model do not change over time. We also write code to fit the model to n segments of a KDE time-series based on a breakpoint analysis. 

(3)The code to analyze the data is contained in the ModelCalibration and RadiocarbonDataCode directories. The analysis was run in "R version 4.2.2 (2022-10-31 ucrt)" with the Rstudio integrated development environment version ``2023.06.0 Build 421." The python code was run using  python version 3.14.6. Python Software Foundation. (2026). Python (Version 3.14.6). https://www.python.org/

(4) Getting started working with the project:

To construct Figs. 4 and 5, construct the KDEs in the RadiocabonDataCode directory. Next, combine the files into the format of the file named NARecessions.csv in the ModelCalibration directory. Next, run the python code in the Model Calibration directory. 

To Construct Figs. 2 and 3, run the Julia code in the ModelCode directory. 
