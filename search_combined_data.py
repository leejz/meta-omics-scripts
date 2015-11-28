#!/usr/bin/env python
"""
--------------------------------------------------------------------------------
Created:  Jackson Lee 1/19/15

This script reads in a combined database data format file and a search file to 
determine functional gene statistics.  

Input file format:
database format:
header	bin_num	function	ontology	organism	taxonomy	Sample_name_RPKM	etc.
scaffold-0_42	1	InsA-like protein	glutathione S-transferase [EC:2.5.1.18]|K00799	Lyngbya sp. PCC 8106	Bacteria;Cyanobacteria;unclassified (derived from Cyanobacteria);Oscillatoriales;unclassified (derived from Oscillatoriales);Lyngbya;Lyngbya sp. PCC 8106;Lyngbya sp. PCC 8106	105.2
...

Query file tab-delimited format:
abbreviation	name	chain or subunit	code	subunit	EC	exclude	nested
DSR	dissimilatory sulfite reductase;sulfite reductase, dissimilatory-type	alpha;beta;A;B	dsr	AB	1.8.99.1		
   
Search pseudocode:
1. Aggregate search columns and extract chain or subunit words
   
2. If KEGG code field is there, set True if matches query
   
3. If EC code field is there, check for matches of the remaining:         
        If no EC field, skip step. 
        If EC mismatch, set fields as False.  
        If EC match found, Set 'maybe'.  
        If 'maybe', check if query chain or subunit exists.
           If either query chain or field exists, search for chain or field match.  
           If either field matches, set True
        If no query field exists, set 'maybe' to True

4. If 3 digit code field is there, check for matches of the remaining:
        If there is no subunit, check only for code
           If code matches, set True
        If subunit field exists, check for subunit and code match
           If subunit and code matches, set True

5. If name field is there, check for name matches of the remaining:
        Set all remainder to False
        If name whole word matches, check if chain exists.
           If chain field matches, set True.
           If chain field mismatch, set False
        If no chain field exists, set True
        
6. If nested field is there, 
        for each nested item, 
           search for named hits:
              set original field to false for all matches
              create a new abbreviation code and field from these an set True if matched
                                                 
7. If anything in exclude column matches in search column, set False
   
8. Set remainder to False
   
NOTES: Careful regex formulation and syntax required for name, subunit, exclude, 
and nested terms terms as these use python's regex search (via Pandas contains 
function). Beware of subset regex hits (e.g. I, III, IIII). Code field is subset 
into [Cc]ode1[subunit1]|[Cc]ode2[subunit1], etc.
      
If a match is found, the data is saved and subdivided by bin.  The ratio is 
calculated with the bin median log2 (RPKM gene / bin RPKM median). This is done 
for all timepoints and the median taken over timepoints. The duplicates are 
averaged.
    
Several tables are output.  The table of bins vs. abbreviations, with averaged 
medians. And a table of gene counts.  This script also outputs a directory of 
search results.    
    
--------------------------------------------------------------------------------
usage:   search_combined_data.py -i data.file -q query.tab -o out.prefix
"""

#-------------------------------------------------------------------------------
#Header - Linkers, Libs, Constants
from string import strip
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import csv
import os
import pandas as pd
import numpy as np	
from re import escape

#-------------------------------------------------------------------------------
#function declarations

#-------------------------------------------------------------------------------
#Body

print "Running..."

if __name__ == '__main__':
    parser = ArgumentParser(usage = "calculate_bin_stats.py -i data.file -o \
out.file",
                            description=__doc__, 
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_filename", action="store", 
                        dest="inputfilename",
                        help="input tab delimited data file")
    parser.add_argument("-q", "--query_filename", action="store", 
                        dest="queryfilename",
                        help="query tab delimited data file")
    parser.add_argument("-o", "--output_filename", action="store", 
                        dest="outputfilename",
                        help="output tab delimited data file")
    options = parser.parse_args()

    mandatories = ["inputfilename", "queryfilename", "outputfilename"]
    
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nError: Missing Arguments\n"
            parser.print_help()
            exit(-1)
    
    print "Reading in and formatting files..."
    inputfilename = options.inputfilename
    queryfilename = options.queryfilename
    outputfilename = options.outputfilename
    outputdirectory = outputfilename.rpartition('.')[0] + '.d'
    
    if not os.path.exists(outputdirectory):
        os.makedirs(outputdirectory)
    
    with open(inputfilename,'U') as infile:
        combined = pd.read_csv(infile, header=0, sep='\t')
        combined = combined.loc[combined.notnull().any(axis=1),:]

    MTnames = combined.columns.tolist()[6:]
    annotnames = ["header", "bin_num", "function", "ontology", "organism", "taxonomy"]
    combined.columns = annotnames + MTnames
    combined.loc[:,"bin_num"] = combined.loc[:,"bin_num"].astype(str)
    
    #bin_list = list(combined.bin_num.unique())
    #bin_list.sort()
    zsum = zip(*[combined.bin_num.groupby(combined.bin_num).count().index, combined.bin_num.groupby(combined.bin_num).count().tolist()])
    zsum.sort(key=lambda x: x[0])
    bin_list = [str(bin) for (bin, count) in zsum if count > 10 and '-x-' not in bin]    
    
    combined.loc[:,'searchcol'] = combined.function.fillna('') +';' + combined.ontology.fillna('')
    combined.loc[:,'EC'] = combined.searchcol.str.extract('\([ ]{0,1}?EC[: ](\d+\.\d+\.\d+\.\d+)[ ]{0,1}?\)')
    comb = combined.searchcol.str.extract('(\w*) [Cc]hain|(\w*) [Ss]ubunit') + ';' + combined.searchcol.str.extract('[Cc]hain (\w*)|[Ss]ubunit (\w*)')
    combined.loc[:,'chain_or_subunit'] = comb.fillna('').sum(axis=1)
        
    with open(queryfilename, 'U') as queryfile:
        #qreader = csv.reader(queryfile, dialect='excel-tab')
        #queries = [for queryline in qreader]
        queries = pd.read_csv(queryfile, header=0, sep='\t', dtype='object')
        queries = queries.loc[queries.notnull().any(axis=1),:]
        
    queries.columns = ["abbreviation", "qname", "chain_or_subunit", "code", "subunit", "EC", "KO", "exclude", "nested"]
    queries.chain_or_subunit = queries.chain_or_subunit.replace(';','|',regex=True)
    #queries.EC = queries.EC.replace(';','|',regex=True)
    queries.EC = queries.EC.replace('\.','\\\.',regex=True)
    #queries.qname = queries.qname.replace(';','|',regex=True)
    #queries.exclude = queries.exclude.replace(';','|',regex=True)
    querynames = queries.abbreviation.tolist()
                            
    print "Searching for functional genes..."    
    
    for qdex in queries.index:
        # 0 abbreviation
        # 1 name
        # 2 chain or subunit
        # 3 code
        # 4 subunit
        # 5 KO
        # 6 EC
        # 7 exclude
        # 8 nested  
        query = queries.iloc[qdex,:]  
        querynotnull = query.notnull()
        queryisnull = query.isnull()
        queryab = query.abbreviation
        combined.loc[:,queryab] = np.nan
        #searching_df = combined[query.abbreviation]
        if querynotnull.KO:
            containsKOpd = combined.loc[combined.ontology.notnull(), 'ontology'].str.contains(query.KO)
            containsKOpd = containsKOpd[containsKOpd == True]
            combined.loc[containsKOpd.index, queryab] =  containsKOpd
        if querynotnull.EC:
            ECsearch_df = combined.loc[combined[queryab].isnull()]
            ECsearch_df.loc[ECsearch_df.EC.notnull(), queryab] = False
            ecmaybe = ECsearch_df.loc[ECsearch_df.EC.str.contains(query.EC) == True]
            if querynotnull.chain_or_subunit or query.notnull().subunit:
                if querynotnull.chain_or_subunit:   
                    qchainquery = query.chain_or_subunit.split('|')
                    ecTF = pd.DataFrame(ecmaybe.chain_or_subunit.str.split(';', expand=True).isin(qchainquery).any(axis=1))                
                if querynotnull.subunit:
                    queryprot = '|'.join([r'\b[' + qcode[0] + qcode[0].upper() + ']' + qcode[1:] + '[' + qsubunit + r']\b' for qcode in query.code.split(';') for qsubunit in query.subunit])
                    if queryisnull.chain_or_subunit:                
                        ecTF = pd.DataFrame(ecmaybe.searchcol.str.contains(queryprot))
                    else:
                        ecTF = ecTF.join(pd.DataFrame(ecmaybe.searchcol.str.contains(queryprot))).any(axis=1)
                combined.loc[ecTF.index, queryab] = ecTF
            else:
                combined.loc[ecmaybe.index, queryab] = True
        if querynotnull.code:
            if queryisnull.subunit:
                querycode = '|'.join([r'\b[' + qcode[0] + qcode[0].upper() + ']' + qcode[1:] + r'\b' for qcode in query.code.split(';')])
                codesearch_df = combined.loc[combined[queryab].isnull(),'searchcol'].str.contains(querycode)
                codesearch_df = codesearch_df[codesearch_df == True]
                combined.loc[codesearch_df.index, queryab] = True
            else:
                queryprot = '|'.join([r'\b[' + qcode[0] + qcode[0].upper() + ']' + qcode[1:] + '[' + qsubunit + r']\b' for qcode in query.code.split(';') for qsubunit in query.subunit])
                codesearch_df = combined.loc[combined[queryab].isnull(),'searchcol'].str.contains(queryprot)
                codesearch_df = codesearch_df[codesearch_df == True]
                combined.loc[codesearch_df.index, queryab] = True		
        if querynotnull.qname:
            namesearch_df = combined.loc[combined[queryab].isnull()]
            combined.loc[namesearch_df.index, queryab] = False
            #ignore case '(?i)'
            queryname = '(?i)' + query.qname
            namematch = namesearch_df[namesearch_df.searchcol.str.contains(queryname).fillna(False)]
            if querynotnull.chain_or_subunit:
                qchainquery = query.chain_or_subunit.split('|')
                namematch = namematch[namematch.chain_or_subunit.str.split(';', expand=True).isin(qchainquery).any(axis=1)]
            combined.loc[namematch.index, queryab] = True	
                                
        if querynotnull.exclude:
            combined.loc[combined.searchcol.str.contains('(?i)' + query.exclude).fillna(False), queryab] = False
                                            
        if querynotnull.nested:
            for nquery in query.nested.split('|'):
                nqueryabbr = queryab + '_' + nquery
                querynames.append(nqueryabbr)
                mnquery = '(?i)' + nquery.replace(';','|')
                allmatches = combined.loc[combined[queryab] == True,'searchcol'].str.contains(mnquery)
                if isinstance(allmatches, pd.DataFrame):
                    allmatches = allmatches.any(axis=1)
                
                combined.loc[allmatches.index, nqueryabbr] = allmatches
                combined.loc[combined[nqueryabbr].isnull(), nqueryabbr] = False
                combined.loc[allmatches[allmatches == True].index, queryab] = False                                
            
        combined.loc[combined[queryab].isnull(), queryab] = False
                                            
#   Search pseudocode:
#   1. Aggregate search columns and extract chain or subunit words
#   
#   2. If KEGG code field is there, set True if matches query
#   
#   3. If EC code field is there, check for matches of the remaining:         
#         If no EC field, skip step. 
#         If EC mismatch, set fields as False.  
#         If EC match found, Set 'maybe'.  
#         If 'maybe', check if query chain or subunit exists.
#            If either query chain or field exists, search for chain or field match.  
#            If either field matches, set True
#         If no query field exists, set 'maybe' to True
#
#   4. If 3 digit code field is there, check for matches of the remaining:
#         If there is no subunit, check only for code
#            If code matches, set True
#         If subunit field exists, check for subunit and code match
#            If subunit and code matches, set True
#
#   5. If name field is there, check for name matches of the remaining:
#         Set all remainder to False
#         If name whole word matches, check if chain exists.
#            If chain field matches, set True.
#            If chain field mismatch, set False
#         If no chain field exists, set True
#   6. If nested field is there, 
#         for each nested item, 
#            search for named hits:
#                set original field to false for all matches
#                create a new abbreviation code and field from these and set True if matched
#                                                 
#   7. If anything in exclude column matches in search column, set False
#   
#   8. Set remainder to False
        
    print "Computing stats and writing..."

    with open(outputfilename, 'w') as outfile:
        writer = csv.writer(outfile, dialect='excel-tab')
        writer.writerow(["bin", 'log 2 median'] + querynames + ['score'] + querynames + ['counts'] + querynames)
        for bin in bin_list:
            #find bins, bin cases: k29.1, 11131.1, k29.13, k29.1-x-11131.1, k29.13-x-11131.1
            #take bins, split on -x-, search for exact match
            #binquery = str(bin).replace('.','\.')
            binquery = str(bin)
            binsearch_df = (combined.bin_num.str.split('-x-', expand=True) == binquery).any(axis=1)
            if binsearch_df.any():
                working_df = combined.loc[binsearch_df, MTnames]
                abbrs_df = combined.loc[binsearch_df, querynames]            
                log_df = np.log(working_df.applymap(float) / working_df.median())
            
                bin_abbr_stats = []
                gene_nums = []
                gene_counts = []
                for abbr in querynames:                   
                    abbr_df = log_df.loc[abbrs_df[abbr]]
                    gene_med = abbr_df.mean()
                    gene_med = gene_med[gene_med != 0].median()
                    bin_abbr_stats.append(gene_med)
                    if np.isnan(gene_med) or np.isinf(gene_med):
                        gene_num = 0
                    elif gene_med < -1:
                        gene_num = 1
                    elif gene_med <= 1:
                        gene_num = 2
                    elif gene_med > 1:
                        gene_num = 3
                    gene_nums.append(gene_num)                    
                    gene_counts.append(len(abbr_df.index))
                    bin_org_df = combined.loc[combined.bin_num == bin, 'organism'].dropna().value_counts()
                    if len(bin_org_df) > 0:
                        bin_org = str(bin_org_df.index[0])
                    else:
                        bin_org = ''
                writer.writerow([bin, bin_org] + bin_abbr_stats + [''] + gene_nums + [''] + gene_counts)     
    
    for abbr in querynames:
        with open(outputdirectory + '/' + abbr + '.search_results.txt', 'w') as abbrfile:
            combined.loc[combined[abbr], annotnames].to_csv(abbrfile, sep="\t", header=True, index=False)
            
    print "Done!"
