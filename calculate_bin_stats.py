#!/usr/bin/env python
"""
--------------------------------------------------------------------------------
Created:  Jackson Lee 12/13/14

This script reads in a database data format file and calculates several statistics.

Input file format:
database format:
header	bin_num	function	ontology	organism	taxonomy	Sample_name_RPKM	etc.
scaffold-0_42	1	InsA-like protein	glutathione S-transferase [EC:2.5.1.18]|K00799	Lyngbya sp. PCC 8106	Bacteria;Cyanobacteria;unclassified (derived from Cyanobacteria);Oscillatoriales;unclassified (derived from Oscillatoriales);Lyngbya;Lyngbya sp. PCC 8106;Lyngbya sp. PCC 8106	105.2
...

--------------------------------------------------------------------------------   
usage:   calculate_bin_stats.py -i data.file -o out.file
"""

#-------------------------------------------------------------------------------
#Header - Linkers, Libs, Constants
from string import strip
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import csv
import pandas as pd
import numpy as np	

#-------------------------------------------------------------------------------
#function declarations

#-------------------------------------------------------------------------------

#Body
print "Running..."

if __name__ == '__main__':
    parser = ArgumentParser(usage = "calculate_bin_stats.py -i data.file -o out.file",
                            description=__doc__, 
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_filename", action="store", 
                        dest="inputfilename",
                        help="input tab delimited data file")
    parser.add_argument("-o", "--output_filename", action="store", 
                        dest="outputfilename",
                        help="output tab delimited data file")
    options = parser.parse_args()

    mandatories = ["inputfilename", "outputfilename"]
    
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nError: Missing Arguments\n"
            parser.print_help()
            exit(-1)

    inputfilename = options.inputfilename
    outputfilename = options.outputfilename
    
    with open(inputfilename,'U') as infile:
        combined = pd.read_csv(infile, header=0, sep='\t')

    combined.columns = ["header", "bin_num", "function", "ontology", "organism", "taxonomy"] + combined.columns.tolist()[6:]
    #bin_list = list(combined.bin_num.unique())
    #bin_list.sort()
    zsum = zip(*[combined.bin_num.groupby(combined.bin_num).count().index, combined.bin_num.groupby(combined.bin_num).count().tolist()])
    zsum.sort(key=lambda x: x[0])
    bin_list = [bin for (bin, count) in zsum if count > 10]    
    
    with open(outputfilename, 'w') as outfile:
        writer = csv.writer(outfile, dialect='excel-tab')
        writer.writerow(["bin"] + combined.columns.tolist()[6:])
        for bin in bin_list:
            working_df = combined[combined.bin_num == bin].iloc[:,6:]
            writer.writerow([str(bin) + ' medians'])
            writer.writerow([bin] + [working_df[colname].median() for colname in working_df.columns.tolist()])
            writer.writerow([str(bin) + ' standard deviation'])
            writer.writerow([bin] + [working_df[colname].std() for colname in working_df.columns.tolist()])
            writer.writerow([str(bin) + ' frequency table'])
            
            minbucket = working_df.min().min() - 0.001
            maxbucket = working_df.median().max() + working_df.std().median()
            buckets = np.arange(minbucket,maxbucket,(maxbucket-minbucket)/101)
            
            frequency_tables = [list(np.histogram(working_df[colname],bins=buckets)[0]) for colname in working_df.columns.tolist()]            
            for i, bucket in enumerate(buckets[:-1]):
                writer.writerow([bucket] + [freqline[i] for freqline in frequency_tables])                        
            
            writer.writerow([str(bin) + ' median normalized frequency table'])
            writer.writerow([str(bin) + ' normalization factor'])
            
            binmed = [working_df.stack().median() / med if med != 0 else 1.0 for med in working_df.median()]
            working_df_norm = working_df * binmed

            writer.writerow([' '] + binmed)
                        
            minbucket = working_df_norm.min().min() - 0.001
            maxbucket = working_df_norm.median().max() + working_df_norm.std().median()
            buckets = np.arange(minbucket,maxbucket,(maxbucket-minbucket)/101)

            frequency_tables = [list(np.histogram(working_df_norm[colname],bins=buckets)[0]) for colname in working_df_norm.columns.tolist()]            
            for i, bucket in enumerate(buckets[:-1]):
                writer.writerow([bucket] + [freqline[i] for freqline in frequency_tables])                        

            
            #ax = working_df.hist()
            #fig = ax.get_figure()
            #fig.savefig('./' + str(bin) + '_histogram.pdf')
                
    print "Done!"
