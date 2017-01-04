#!/usr/bin/env python
"""
--------------------------------------------------------------------------------
Created:  Jackson Lee 9/27/14

This script reads in a tab delimited file of annotations and querys the KEGG 
REST API to parse back the original KO ortholog of the entry.  If the line 
contains an EC reference, the script will first query each line in KEGG REST as:
http://rest.kegg.jp/find/genes/each+query+term+and+EC+d.d.d.d

e.g. 1-deoxy-D-xylulose 5-phosphate synthase (EC 2.2.1.7)
http://rest.kegg.jp/find/genes/1-deoxy-D-xylulose+5-phosphate+synthase+2.2.1.7

and save the output in query order with the format:
gmx:100301901	DXS1; 1-deoxy-D-xylulose 5-phosphate synthase 1; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
pvu:PHAVU_009G095900g	hypothetical protein; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
pvu:PHAVU_006G159900g	hypothetical protein; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
pvu:PHAVU_003G099600g	hypothetical protein; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
pvu:PHAVU_003G148900g	hypothetical protein; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
pvu:PHAVU_003G287800g	hypothetical protein; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
pvu:PHAVU_003G287900g	hypothetical protein; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
mtr:MTR_2g020590	1-deoxy-D-xylulose 5-phosphate synthase; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
mtr:MTR_3g107740	hypothetical protein; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
mtr:MTR_4g118640	1-deoxy-D-xylulose 5-phosphate synthase; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
mtr:MTR_8g068300	1-Deoxy-D-xylulose 5-phosphate synthase; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
mtr:MTR_8g068270	1-Deoxy-D-xylulose 5-phosphate synthase; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]
mtr:MTR_8g068280	1-Deoxy-D-xylulose 5-phosphate synthase; K01662 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]

This output will then be read in and queried for the exact search term and grab by regex
any Kterms in the string.  These terms are aggregated and the top hit and score written out
regex: "1-deoxy-D-xylulose 5-phosphate synthase;" and "EC:2.2.1.7" for "; K\d{5}" 

result:KO1662
   
   Input file format:
(3R)-hydroxymyristoyl-[ACP] dehydratase (EC 4.2.1.-)
(R)-citramalate synthase (EC 2.3.1.182)
(S)-2-haloacid dehalogenase I (EC 3.8.1.2)
(S)-22C3-di-O-geranylgeranylglyceryl phosphate synthase
(S)-3-O-geranylgeranylglyceryl phosphate synthase
(Y14336) putative extracellular protein containing predicted 35aa signal peptide
1-acyl-sn-glycerol-3-phosphate acyltransferase (EC 2.3.1.51)
1-aminocyclopropane-1-carboxylate deaminase (EC 3.5.99.7)
1-deoxy-D-xylulose 5-phosphate reductoisomerase (EC 1.1.1.267)
1-deoxy-D-xylulose 5-phosphate synthase (EC 2.2.1.7)
   
   Output
   a translation table of terms and the KEGG REST output
   
   1-deoxy-D-xylulose 5-phosphate synthase (EC 2.2.1.7)\tK01662\t5\t5

--------------------------------------------------------------------------------     
usage:   query_ko_from_rast_annotations.py -i in.file -d out.directory -o output.file
"""

#-------------------------------------------------------------------------------
#
#http thread pool code from: http://stackoverflow.com/questions/2632520/what-is-the-fastest-way-to-send-100-000-http-requests-in-python 
#-------------------------------------------------------------------------------
#Header - Linkers, Libs, Constants
from string import strip
import os
import re
import collections
from argparse import ArgumentParser, RawDescriptionHelpFormatter
#import requests

from urlparse import urlparse
from threading import Thread
import httplib
import sys
from Queue import Queue

#-------------------------------------------------------------------------------
#function declarations

def doWork():
    while not exitapp:
        id, urlstring, queryline = q.get()
        url = urlparse(urlstring)
        if id % 100 == 0:
            print 'Query: HTTP Thread: ' + str(id) + ' started.'
        try:
            conn = httplib.HTTPConnection(url.netloc)
            conn.request("GET", url.path)
            res = conn.getresponse()

            if res.status == 200:
                with open(outputdirectory + '/' + str(id) + '.KEGG_REST.txt', 'w') as restfile:
                    restfile.write(res.read())
                restfile.close()
                #print 'Thread: ' + str(id) + ' Query: ' + urlstring + ' ..... ' + res.reason + '\n'
                searchfile.write(str(id) + '\t' + queryline + '\n')
            else:
                print 'HTTP error, Thread: ' + str(id) + ' with error: ' + res.reason
                logfile.write(str(id) + '\t' + urlstring + '\t' + res.reason + '\n')
                raise
        except:
            print 'Thread: ' + str(id) + '. Error. ' 
            print sys.exc_info()[0]
            
        q.task_done()

#-------------------------------------------------------------------------------
#Body
print "Running..."

if __name__ == '__main__':
    parser = ArgumentParser(usage = "query_ko_from_rast_annotations.py -i \
in.file -d out.directory -o output.file",
                            description=__doc__, 
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_file", action="store", 
                        dest="inputfilename", help="text input file")
    parser.add_argument("-o", "--output_filename", action="store", 
                        dest="outputfilename", help="text output file")
    parser.add_argument("-d", "--output_directory", action="store", 
                        dest="outputdirectory", help="text output file")
    options = parser.parse_args()

    mandatories = ["outputfilename","outputdirectory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nError: Missing Arguments\n"
            parser.print_help()
            exit(-1)
    outputdirectory = options.outputdirectory            
    
    ec_regex = '\d\.[\d\-]*\.[\d\-]*\.[\d\-]*'
            
    #allow for only querying the KEGG REST API once
    if options.__dict__['inputfilename']:
        if not os.path.exists(outputdirectory):
            os.makedirs(outputdirectory)
        else:
            print "\nError: Directory exists!\n"
            parser.print_help()
            exit(-1)

        print "Querying KEGG REST API Service..."     
        inputfilename = options.inputfilename    

        api_template = "http://rest.kegg.jp/find/genes/"
        infile_list = []
        with open(inputfilename,'U') as infile:
            infile_list = [line.strip() for line in infile]
        infile.close()
        # replace all 2C, %3B
        infile_list = [line.replace('2C',',') for line in infile_list]
        infile_list = [line.replace('%3B',';') for line in infile_list]

        urlpool = []
        for line in infile_list:
            if re.search('EC ' + ec_regex, line) != None:
                #format string for search
                query = line.strip()
                #remove ec
                ecnum_list = re.findall('EC ' + ec_regex,query)
                ecnum_list = [ecline[3:] for ecline in ecnum_list]
                query = re.sub(' \(EC ' + ec_regex + '\)', '', query)
                #remove url syntax issues
                query = query.replace('+','')
                query = query.replace('/',' ')
                query = query.replace('@',' ')
                query = query.replace(';',' ')
                query = query.replace ('  ',' ')
                #query += ' '.join(ecnum_list)

                #urlstring = api_template + query

                #form url, query, and write file
                querylist = filter(None, query.split(' ') + ecnum_list)                                        
                urlstring = api_template + '+'.join(querylist)
                #catch case of '+-' '+)+' '+(+' and convert url encoding
                urlstring = urlstring.replace('+-','+')
                #urlstring = urlstring.replace('+)+','+')
                #urlstring = urlstring.replace('+(+','+')
                #urlstring = urlstring.replace(' ','%20')                                
                urlstring = urlstring.replace('(','%28')                                
                urlstring = urlstring.replace(')','%29')
                urlpool.append([urlstring, line])    
                    
                #    print 'Query: ' + urlstring
                #    r = requests.get(urlstring)
                #    if r.raise_for_status() == None:
                #        with open(outputdirectory + '/' + str(i) + '.KEGG_REST.txt', 'w') as restfile:
                #            restfile.write(r.text)
                #        restfile.close()
                #        searchfile.write(str(i) + '\t' + line + '\n')
                #    else:
                #        print 'Response error raised.  Exiting'
                #        exit(-1)
 
        #setup threading for http requests and run connections
        concurrent = 100
        exitapp = False
        
        with open(outputdirectory + '/searchlist.txt', 'w') as searchfile, open(outputdirectory + '/errorlog.txt','w') as logfile:           
            q = Queue(concurrent * 2)
            for i in range(concurrent):
                t = Thread(target=doWork)
                t.daemon = True
                t.start()
            try:
                for id, urlentry in enumerate(urlpool):
                    q.put([id] + urlentry)
                q.join()
            except KeyboardInterrupt:
                exitapp = True
                sys.exit(1)                            
                logfile.close()
        searchfile.close()
        logfile.close()
        
    
    print "Parsing REST files and writing..."
            
    outputfilename = options.outputfilename
    outfile = open(outputdirectory + '/' + outputfilename, 'w')
    with open(outputdirectory + '/searchlist.txt','U') as searchfile:
        for line in searchfile:
            i, query = line.strip().split('\t')
            #form string for search
            ecnum_list = re.findall('EC ' + ec_regex,query)
            ecnum_list = [ecline[3:] for ecline in ecnum_list]
            #querystring = '\t' + re.sub(' \(EC ' + ec_regex + '\)', '', query)
            querystrings = [re.sub(' \(EC ' + ec_regex + '\)', '', querystring).lower() for querystring in re.split(' / | @ |; ', query)]
            ecstring = '(EC:' + ' '.join(ecnum_list) + ');'
            ko = []            
            
            with open(outputdirectory + '/' + str(i) + '.KEGG_REST.txt', 'U') as inrestfile:
                for restline in inrestfile:
                    restline = restline.strip()
                    #if querystring == 'Chlorophyllide reductase subunit BchZ':
                    #    print querystring, ecstring
                    #    print querystring in restline, ecstring in restline
                    # if the enzyme search string and the modified ec string and a KEGG KO number are in the rest output, record the KO number
                    if all(querystring in restline.lower() for querystring in querystrings) and all(ecterm in restline for ecterm in ecnum_list) and re.search(r'; K\d{5}', restline) != None:
                        ko.append(re.search(r'; K\d{5}', restline).group(0)[2:])
            inrestfile.close()
            #determine and record the most common KO number and how common it was
            counter= collections.Counter(ko)
            if len(counter) > 0:    
                outfile.write(query + '\t' + counter.most_common(1)[0][0] + '\t' + str(counter.most_common(1)[0][1]) + '\t' + str(sum(counter.values())) + '\t' + querystring + '\t' + ecstring + '\n')
            else:
                outfile.write(query + '\t\n')
    outfile.close()

    print "Done!"
