You have been provided with raw data extracted from the experiment 
and a series of Matlab functions and scripts to process the raw data. 
This processing software provides the uncorrected data as an output. 
In order to correct the data you need to include the balance raw data
files together with the zero measurement data. Please follow the steps 
below.

----------------------------------------------------------------
---------------------------STEPS--------------------------------
----------------------------------------------------------------

1. Open 'main_AE4115labExercise_balance.m'

2. Modify the diskPath (line 13).

3. Enter the filename(s) of the 'raw_..' balance data 
(lines 26 and 30) separated with commas. Don't forget
to include the extension(.txt). Don't use the 'unc_..' and 'cor_..' 
files, only the 'raw_..' files. The 'unc_..' file includes live 
processed data from the wind-tunnel computer. The results from the 
matlab script will be more accurate (mostly due to the way the zero 
measurements are treated). The 'cor_..' file includes live corrected 
processed data from the wind-tunnel computer. These corrections are 
not correct nor complete because some corrections are missing and not 
all inputs provided to the wind-tunnel computer have been set up. 
You will implement your own code to correct the data that are the 
output of the Matlab script.

4. Enter the filename(s) of the zero measurements (line 39). 
You have to enter one file for each balance data file that you have. 
If you only have one zero measurement file, you have to
write the name of the file X times.

5. There is a function that writes the results in an 'xls'
file just in case that you prefer to work with Excel. (Ignore the 
empty sheets sheet1, sheet2, sheet3)
Use of the excel sheet is optional, you can work with the output 
data in Matlab. Just write your code after line 73.

6. Run the script. The output data will be created as a
structure type with the name 'BAL' and 'PRS'. 

I strongly recommend you to have a look at all the functions
in order to understand the output number that is being
provided to you.