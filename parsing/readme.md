INSTRUCTIONS FOR RUNNING GO SHORTLISTING PIPELINE TOOL

===================================================================================

Step 1: Prepare GO list for analysis using PANTHER database
Go to http://pantherdb.org/tools/compareToRefList.jsp. Enter the list of differentially expressed (DE) genes you wish to analyze as the analyzed list and the list of all genes on the microarray chip as the reference list. Both should be text files with genes in a single column to be accepted by the database. Select 'GO biological process complete' for the annotation dataset and run the overrepresentation analysis. Download the resulting JSON file and save as species_timept.json, for the respective species the genes are from and time point that the DE genes were measured in the experiment, in the appropriate location. 

Step 2: Run representation filtering in Python
On the command line, in the appropriate directory, run the following:
python3 shortlisting.py species timept num_in_list rep
where num_in_list=the number of genes in the analyzed list.
The command line will prompt you to enter biases on representation metrics Ra and Rn. These should be numbers between 0 and 1 as prescribed by the analysis you wish to perform-- see paper.

Step 3: Prepare similarity scores using NaviGO database
Open 'species_timept_filtered.xlsx'. Go to https://kiharalab.org/web/navigo/views/goset.php. Copy the list of GO terms from the Excel sheet and enter them into the database under Input GO Terms and submit. The database will produce a list of similarity scores which is downloadable as a CSV file. Download the file and save as an Excel file using the name 'species_timept_sim.xlsx'.

Step 4: Run similarity filtering in Python
On the command line, in the appropriate directory, run the following:
python3 shortlisting.py species timept num_in_list sim

Step 5: View completed shortlist by opening Excel file 'species_timept_shortlist.xlsx'.

===================================================================================

Please reach out to esphilli@ucsc.edu with any questions. 