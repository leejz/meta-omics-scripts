#!/usr/bin/env python
"""
--------------------------------------------------------------------------------
Created: Jackson Lee 9/29/14

This script reads in a tab delimited combined coverage file from 
consolidate_coverage.py of phylogeny, protein classifications, and annotations 
uses the biopython KGML libraries to generate graphics of KEGG pathways.  

Input file:
Phylogeny	Organism	Protein Classification	RPKM1	RPKM2 ...
2	H. Monster	Function|K00003	4211.629513	...
2	H. Monster	Function|K00012	2752.574388
3	...		...
   
   Output
A series of mapping files for each bin over each time point

--------------------------------------------------------------------------------
usage:   generate_kgml_graphics.py -i in.file -d out.directory

"""

#-------------------------------------------------------------------------------
#
#
#Code from: http://armchairbiology.blogspot.co.uk/2013/02/keggwatch-part-iii.html
#Generating KEGG maps example 2
#
#import KGML_parser
#from KGML_scrape import retrieve_KEGG_pathway
#from KGML_vis import KGMLCanvas
#from Bio.Graphics.ColorSpiral import ColorSpiral
 
## Get the ko03070 map from KEGG, and write it out to file, visualised as
## the .png, and as the elements from the KGML file
#pathway = retrieve_KEGG_pathway('ko03070')
#kgml_map = KGMLCanvas(pathway, show_maps=True)
# 
## Let's use some arbitrary colours for the orthologs
#cs = ColorSpiral(a=2, b=0.2, v_init=0.85, v_final=0.5,
#jitter=0.03)
## Loop over the orthologs in the pathway, and change the
## background colour
#orthologs = [e for e in pathway.orthologs]
#for o, c in zip(orthologs,
#cs.get_colors(len(orthologs))):
#for g in o.graphics:
#g.bgcolor = c
## Default settings are for the KGML elements only
#kgml_map.draw('ex2_kgml_render.pdf')
# 
## We need to use the image map, and turn off the KGML elements, to see
## only the .png base map. We could have set these values on canvas
## instantiation
#kgml_map.import_imagemap = True
#kgml_map.show_maps = False
#kgml_map.show_orthologs = False
#kgml_map.draw_relations = False
#kgml_map.show_compounds = False
#kgml_map.show_genes = False
#kgml_map.draw('ex2_png_render.pdf')
# 
## And rendering elements as an overlay
#kgml_map.show_compounds = True
#kgml_map.show_genes = True
#kgml_map.show_orthologs = True
#kgml_map.draw('ex2_overlay_render.pdf')
#
#Generating KEGG Maps example 3
#
#import KGML_parser
#from KGML_scrape import retrieve_KEGG_pathway
#from KGML_vis import KGMLCanvas
# 
## Get list of pathway elements to enhance
#glyc_path = retrieve_KEGG_pathway('ko00010')
#tca_path = retrieve_KEGG_pathway('ko00020')
#enhance_list = []
#for pathway in (glyc_path, tca_path):
#for e in pathway.entries.values():
#enhance_list.extend(e.name.split())
#enhance_list = set(enhance_list)
# 
## Get the pathway we want to render, and make all the lines
## that are also in glycolysis or TCA pathways thicker
#met_pathway = retrieve_KEGG_pathway('ko01100')
#mod_list = [e for e in met_pathway.entries.values() if \
#len(set(e.name.split()).intersection(enhance_list))]
#for e in mod_list:
#for g in e.graphics:
#g.width = 10
#kgml_map = KGMLCanvas(met_pathway, show_maps=True)
#kgml_map.draw('ex3_thick.pdf')
## Thin out any lines that aren't in the glycolysis/TCA pathways
#mod_list = [e for e in met_pathway.entries.values() if \
#not len(set(e.name.split()).intersection(enhance_list)) \
#and e.type != 'map']
#for e in mod_list:
#for g in e.graphics:
#g.width = .4
#kgml_map.draw('ex3_thin.pdf')
# 
## Or turn them grey, maybe:
#for e in mod_list:
#for g in e.graphics:
#g.fgcolor = '#CCCCCC'
#kgml_map.draw('ex3_grey.pdf')
#
#Code from: http://stackoverflow.com/questions/452074/creating-a-gradient-fill-in-a-pdf-file-using-reportlab
#
#from reportlab.pdfgen.canvas import Canvas
#from reportlab.lib.colors import red, yellow, green
#from reportlab.lib.units import mm
#
#c = Canvas("gradient.pdf")
#
## Linear gradient with the endpoints extending over the page.
#c.linearGradient(105*mm, 200*mm, 180*mm, 100*mm, (red, yellow))
#c.drawString(5*mm, 290*mm, "c.linearGradient(105*mm, 200*mm, 180*mm, 100*mm, (red, yellow))")
#c.line(105*mm, 200*mm, 180*mm, 100*mm)
#c.showPage()
#
## Linear gradient constrained within the endpoints.
#c.linearGradient(105*mm, 200*mm, 180*mm, 100*mm, (red, yellow), extend=False)
#c.drawString(5*mm, 290*mm, "c.linearGradient(105*mm, 200*mm, 180*mm, 100*mm, (red, yellow), extend=False)")
#c.line(105*mm, 200*mm, 180*mm, 100*mm)
#c.showPage()
#
## Linear gradient with multiple stops.
#c.linearGradient(105*mm, 200*mm, 180*mm, 100*mm, (red, yellow, green), (0, 0.8, 1), extend=False)
#c.drawString(5*mm, 290*mm, "c.linearGradient(105*mm, 200*mm, 180*mm, 100*mm, (red, yellow, green), (0, 0.8, 1), extend=False)")
#c.line(105*mm, 200*mm, 180*mm, 100*mm)
#c.line(141*mm, 102*mm, 189*mm, 138*mm)
#c.showPage()
#
## Radial gradient with the endpoint extending over the page.
#c.radialGradient(105*mm, 200*mm, 60*mm, (red, yellow))
#c.drawString(5*mm, 290*mm, "c.radialGradient(105*mm, 200*mm, 60*mm, (red, yellow))")
#c.circle(105*mm, 200*mm, 60*mm)
#c.showPage()
#
## Radial gradient constrained within the circle.
#c.radialGradient(105*mm, 200*mm, 60*mm, (red, yellow), extend=False)
#c.drawString(5*mm, 290*mm, "c.radialGradient(105*mm, 200*mm, 60*mm, (red, yellow), extend=False)")
#c.circle(105*mm, 200*mm, 60*mm)
#c.showPage()
#
## Radial gradient with multiple stops.
#c.radialGradient(105*mm, 200*mm, 60*mm, (red, yellow, green), (0, 0.8, 1))
#c.drawString(5*mm, 290*mm, "c.radialGradient(105*mm, 200*mm, 60*mm, (red, yellow, green), (0, 0.8, 1))")
#c.circle(105*mm, 200*mm, 48*mm)
#c.circle(105*mm, 200*mm, 60*mm)
#c.showPage()
#
#c.save()
#
#-------------------------------------------------------------------------------

#Header - Linkers, Libs, Constants
from string import strip
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import pandas as pd
import os
import bisect
from numpy import log10, arange

import KGML_parser
from KGML_scrape import retrieve_KEGG_pathway
from KGML_vis import KGMLCanvas
from Bio.Graphics import ColorSpiral

# List of 2010 IDs for metabolic pathways
metabolic = ["ko00010", "ko00020", "ko00030", "ko00040", "ko00051", "ko00052",
             "ko00053", "ko00061", "ko00062", "ko00071", "ko00072", "ko00100",
             "ko00120", "ko00121", "ko00130", "ko00140", "ko00190", "ko00195",
             "ko00196", "ko00230", "ko00231", "ko00232", "ko00240", "ko00250",
             "ko00253", "ko00260", "ko00270", "ko00280", "ko00281", "ko00290",
             "ko00300", "ko00310", "ko00311", "ko00312", "ko00330", "ko00331",
             "ko00340", "ko00350", "ko00351", "ko00360", "ko00361", "ko00362",
             "ko00363", "ko00364", "ko00380", "ko00400", "ko00401", "ko00402",
             "ko00410", "ko00430", "ko00440", "ko00450", "ko00460", "ko00471",
             "ko00472", "ko00473", "ko00480", "ko00500", "ko00510", "ko00511",
             "ko00512", "ko00513", "ko00514", "ko00520", "ko00521", "ko00522",
             "ko00523", "ko00524", "ko00531", "ko00532", "ko00533", "ko00534",
             "ko00540", "ko00550", "ko00561", "ko00562", "ko00563", "ko00564",
             "ko00565", "ko00590", "ko00591", "ko00592", "ko00600", "ko00601",
             "ko00603", "ko00604", "ko00620", "ko00621", "ko00622", "ko00623",
             "ko00624", "ko00625", "ko00626", "ko00627", "ko00630", "ko00633",
             "ko00640", "ko00642", "ko00643", "ko00650", "ko00660", "ko00670",
             "ko00680", "ko00710", "ko00720", "ko00730", "ko00740", "ko00750",
             "ko00760", "ko00770", "ko00780", "ko00785", "ko00790", "ko00791",
             "ko00830", "ko00860", "ko00900", "ko00901", "ko00902", "ko00903",
             "ko00904", "ko00905", "ko00906", "ko00908", "ko00909", "ko00910",
             "ko00920", "ko00930", "ko00940", "ko00941", "ko00942", "ko00943",
             "ko00944", "ko00945", "ko00950", "ko00960", "ko00965", "ko00966",
             "ko00970", "ko00980", "ko00981", "ko00982", "ko00983", "ko01040",
             "ko01051", "ko01053", "ko01055", "ko01056", "ko01057", "ko01058",
             "ko01100", "ko01110", "ko01120", "ko04070"]

# List of 2010 IDs for non-metabolic pathways
non_metabolic = ["ko02010", "ko02020", "ko02030", "ko02040", "ko02060",
                 "ko03008", "ko03010", "ko03013", "ko03015", "ko03018",
                 "ko03020", "ko03022", "ko03030", "ko03040", "ko03050",
                 "ko03060", "ko03070", "ko03320", "ko03410", "ko03420",
                 "ko03430", "ko03440", "ko03450", "ko04010", "ko04011",
                 "ko04012", "ko04013", "ko04020", "ko04060", "ko04062",
                 "ko04070", "ko04075", "ko04080", "ko04110", "ko04111",
                 "ko04112", "ko04113", "ko04114", "ko04115", "ko04120",
                 "ko04122", "ko04130", "ko04140", "ko04141", "ko04142",
                 "ko04144", "ko04145", "ko04146", "ko04150", "ko04210",
                 "ko04260", "ko04270", "ko04310", "ko04320", "ko04330",
                 "ko04340", "ko04350", "ko04360", "ko04370", "ko04380",
                 "ko04510", "ko04512", "ko04514", "ko04520", "ko04530",
                 "ko04540", "ko04610", "ko04612", "ko04614", "ko04620",
                 "ko04621", "ko04622", "ko04623", "ko04626", "ko04630",
                 "ko04640", "ko04650", "ko04660", "ko04662", "ko04664",
                 "ko04666", "ko04670", "ko04672", "ko04710", "ko04711",
                 "ko04712", "ko04720", "ko04722", "ko04730", "ko04740",
                 "ko04742", "ko04744", "ko04745", "ko04810", "ko04910",
                 "ko04912", "ko04914", "ko04916", "ko04920", "ko04930",
                 "ko04940", "ko04950", "ko04960", "ko04961", "ko04962",
                 "ko04964", "ko04966", "ko04970", "ko04971", "ko04972",
                 "ko04973", "ko04974", "ko04975", "ko04976", "ko04977",
                 "ko04978", "ko05010", "ko05012", "ko05014", "ko05016",
                 "ko05020", "ko05100", "ko05110", "ko05111", "ko05120",
                 "ko05130", "ko05131", "ko05140", "ko05142", "ko05143",
                 "ko05144", "ko05145", "ko05146", "ko05150", "ko05152",
                 "ko05160", "ko05162", "ko05200", "ko05210", "ko05211",
                 "ko05212", "ko05213", "ko05214", "ko05215", "ko05216",
                 "ko05217", "ko05218", "ko05219", "ko05220", "ko05221",
                 "ko05222", "ko05223", "ko05310", "ko05320", "ko05322",
                 "ko05323", "ko05330", "ko05332", "ko05340", "ko05410",
                 "ko05412", "ko05414", "ko05416"]

#all_kegg = metabolic + non_metabolic
#essential
all_kegg = ["ko00010", "ko00020", "ko00030", "ko00040", "ko00051", "ko00052", "ko00053", "ko00061", "ko00071", "ko00190", "ko00195", "ko00196", "ko00230", "ko00240", "ko00250", "ko00260", "ko00270", "ko00500", "ko00510", "ko00520", "ko00562", "ko00620", "ko00625", "ko00630", "ko00640", "ko00650", "ko00660", "ko00680", "ko00710", "ko00720", "ko00910", "ko00920", "ko01100", "ko01110", "ko01120", "ko02010", "ko02020", "ko02060", "ko03070", "ko04710"]

#for bin 27
#all_kegg = ["ko00010", "ko00020", "ko00030", "ko00190", "ko00195", "ko00620", "ko00630", "ko00640", "ko00650", "ko00660", "ko00680", "ko00720", "ko00910", "ko00920", "ko01100", "ko01110", "ko01120", "ko02010", "ko02020", "ko03070", "ko04122"]

#bare set
#all_kegg = ["ko00010", "ko00020", "ko01100", "ko01110", "ko01120"]

#-------------------------------------------------------------------------------

#function declarations


             

#-------------------------------------------------------------------------------
#Body
print "Running..."

if __name__ == '__main__':
    parser = ArgumentParser(usage = "generate_kgml_graphics.py -i in.file -d \
out.directory",                  
                  description=__doc__, 
                  formatter_class=RawDescriptionHelpFormatter)    
    parser.add_argument("-i", "--input_file", action="store", 
                  dest="inputfilename", help="text input file")
    parser.add_argument("-d", "--output_directory", action="store", 
                  dest="outputdirectory", help="text output file")
    parser.add_argument("-K", "--KEGG_directory", action="store", 
                  dest="KEGGdirectory", help="path to KEGG kgml files")                  
    options = parser.parse_args()

    mandatories = ["inputfilename","outputdirectory", "KEGGdirectory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nError: Missing Arguments\n"
            parser.print_help()
            exit(-1)
    outputdirectory = options.outputdirectory            
    inputfilename = options.inputfilename        
    keggdir = options.KEGGdirectory
            
    if not os.path.exists(outputdirectory):
        os.makedirs(outputdirectory)
    else:
        print "\nError: Directory exists!\n"
        parser.print_help()
        exit(-1)

    print "Reading in datafile..."
    
    with open(inputfilename,'U') as infile:
        combined = pd.read_csv(infile, header=0, sep='\t')
    infile.close()
    combined.columns = ["Phylogeny", "Organism", "Protein Classification"] + combined.columns.tolist()[3:]
    combined["Protein Classification"] = combined["Protein Classification"].str.replace('^.*\|', '')
    rpkm_columns = combined.columns[3:]
    log10_columns = [column_name + '_log10' for column_name in rpkm_columns]
    combined[log10_columns] = combined[rpkm_columns].applymap(lambda x: log10(float(x)) if float(x) > 0 else 0)
    bin_list = list(combined.Phylogeny.unique())
    bin_list.sort()
    
    #cs = ColorSpiral(a=2, b=0.2, v_init=0.85, v_final=0.5, jitter=0.03)        
    print "Generating graphics..."     
    for bin in bin_list:
        working_df = combined[combined.Phylogeny == bin]
        os.makedirs(outputdirectory + '/' + str(bin))
        #.reindex(index='Protein Classification') 
        for timepoint, label in zip(log10_columns,rpkm_columns):
            # find rpkm ranges and set color palette
            min_rpkm = working_df[working_df[timepoint] != 0][timepoint].min()
            max_rpkm = working_df[working_df[timepoint] != 0][timepoint].max() 
            cutoff_rpkm = working_df[working_df[timepoint] != 0][timepoint].median()
            color_range = arange(min_rpkm, max_rpkm, (max_rpkm-min_rpkm)/100)         
            color_dict = ColorSpiral.get_color_dict(color_range, a=6, b=0.7, v_init=0.7, v_final=0.55, jitter=0.00)                          
            print 'Generating ' + outputdirectory + '/' + str(bin) + '/' + str(bin) + '.' + label
            for map in all_kegg:                           
                outfilename = outputdirectory + '/' + str(bin) + '/' + str(bin) + '.' + label + '.' + map + '.pdf'
                #print 'Opening ' + keggdir + '/' + map + '.kgml'
                pathway = KGML_parser.read(open(keggdir + '/' + map + '.kgml', 'U'))
                kgml_map = KGMLCanvas(pathway, show_maps=False)
                kgml_map.fontsize = 9
                special_maps = ['ko01100','ko01110','ko01120']
                
                if pathway.name.split('path:')[1] in special_maps:
                    entries = [e for e in pathway.orthologs]
                    for entry in entries:
                        ename = entry.name.split('ko:')[1:]
                        ename = [i[:6].lower() for i in ename]
                        erpkm = working_df.loc[working_df["Protein Classification"].isin(ename),label].sum()
                        if erpkm >= 0:
                            erpkm = log10(erpkm)
                        
                        if erpkm < min_rpkm:
                            #print oname
                            for g in entry.graphics:
                                g.fgcolor = '#CCCCCC'
                                g.width = .4
                        else:
                            for g in entry.graphics:
                                g.width = 2
                        if erpkm > cutoff_rpkm:
                            for g in entry.graphics:
                                g.width = 10
                        
                    kgml_map.show_colorbar_legend = False
                    kgml_map.import_imagemap = False
                    kgml_map.show_maps = True
                    
                    
                else:
                    kgml_map.set_colorbar_legend(minmax=['1e%.2f' % min_rpkm,'1e%.2f' % max_rpkm], wh_dims = [60.0, 5.0], xypos= [35.0, 5.0], color_dict=color_dict)
                    orthologs = [e for e in pathway.orthologs]
                    for ortholog in orthologs:
                        oname = ortholog.name.split('ko:')[1:]
                        oname = [i[:6].lower() for i in oname]
                        orpkm = working_df.loc[working_df["Protein Classification"].isin(oname),label].sum()
                        if orpkm != 0:
                            orpkm = log10(orpkm)
                        if orpkm > max_rpkm:
                            orpkm = max_rpkm
                        if orpkm <= 0:
                            orpkm = min_rpkm
                        if bisect.bisect_left(color_range, orpkm) > len(color_range)-1:
                            ocolor = color_dict[color_range[-1]]
                        else:
                            ocolor = color_dict[color_range[bisect.bisect_left(color_range, orpkm)]]
                        for element in ortholog.graphics:
                            element.bgcolor = ocolor
                    
                    # set figure display attributes
                    kgml_map.import_imagemap = True
                    #kgml_map.show_maps = True
                    kgml_map.show_maps = False
                    #kgml_map.show_orthologs = False
                    kgml_map.draw_relations = False
                    kgml_map.show_compounds = False
                    #kgml_map.show_genes = False
 
                    # And rendering elements as an overlay
                    #kgml_map.show_compounds = True
                    kgml_map.show_genes = True
                    kgml_map.show_orthologs = True

                # Default settings are for the KGML elements only
                kgml_map.draw(outfilename)           

    print "Done!"
