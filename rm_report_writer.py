#!/usr/bin/python36

##################################################################
# Import modules                                                 #
##################################################################
import re
import glob
import sys

##################################################################
# Define global varialbes                                        #
##################################################################
indices = (["sines", "alus", "b2b4","ids", "mirs", "lines", "line1", "line2", "l3cr1","ltr",
            "ervl","ervl_malrs", "erv_classI", "erv_classII", "DNA", "hat_charlie", "tcmar", "unclassified",
            "small_rna", "satellites", "simple","low_complex" ])


##################################################################
# Program definitions                                            #
##################################################################

def process_file(file_name):
    num_elements = {}
    len_occupied = {}

    sequences = 0
    total_length = 0
    total_length_excl_NX_runs = 0
    bases_masked = 0
    interspersed_repeats = 0

    def populate_dictionaries(key, value1, value2):
        if not value2.isdigit():
            value2 = value2[:-2]
        num_elements[key] = int(value1)
        len_occupied[key] = int(value2)
        return

    f = open(file_name, "r")
    for line in f:
        line.strip("\n")
        line = re.sub(" +", " ", line)
        words = line.split(" ")
        if words[0] == "":
            del words[0]
        if words[0] == "sequences:":
            sequences = int(words[-1])
        elif words[0:2] == ["total", "length:"]:
            total_length = int (words[2])
            total_length_excl_NX_runs =  int(words[4][1:])
        elif words[0:2] == ["bases", "masked:"]:
            bases_masked = int(words[2])
        elif words[0] == "SINEs:":
            populate_dictionaries("sines", words[1], words[2])
        elif words[0] == "Alu/B1":
            populate_dictionaries("alus", words[1], words[2])
        elif words[0] == "B2-B4":
            populate_dictionaries("b2b4", words[1], words[2])
        elif words[0] == "IDs":
            populate_dictionaries("ids", words[1], words[2])
        elif words[0] == "MIRs":
            populate_dictionaries("mirs", words[1], words[2])
        elif words[0] == "LINEs:":
            populate_dictionaries("lines", words[1], words[2])
        elif words[0] == "LINE1":
            populate_dictionaries("line1", words[1], words[2])
        elif words[0] == "LINE2":
            populate_dictionaries("line2", words[1], words[2])
        elif words[0] == r"L3/CR1":
            populate_dictionaries("l3cr1", words[1], words[2])
        elif words[0:2] == ["LTR","elements:"]:
            populate_dictionaries("ltr", words[2], words[3])
        elif words[0] == "ERVL":
            populate_dictionaries("ervl", words[1], words[2])
        elif words[0] == "ERVL-MaLRs":
            populate_dictionaries("ervl_malrs", words[1], words[2])
        elif words[0] == "ERV_classI":
            populate_dictionaries("erv_classI", words[1], words[2])
        elif words[0] == "ERV_classII":
            populate_dictionaries("erv_classII", words[1], words[2])
        elif words[0:2] == ["DNA","elements:"]:
            populate_dictionaries("DNA", words[2], words[3])
        elif words[0] == "hAT-Charlie":
            populate_dictionaries("hat_charlie", words[1], words[2])
            num_elements["hat_charlie"] = int(words[1])
            len_occupied["hat_charlie"] = int (words[2])
        elif words[0] == "TcMar-Tigger":
            populate_dictionaries("tcmar", words[1], words[2])
        elif words[0] == "Unclassified:":
            populate_dictionaries("unclassified", words[1], words[2])
        elif words[0:2] == ["Small", "RNA:"]:
            populate_dictionaries("small_rna", words[2], words[3])
        elif words[0] == "Satellites:":
            populate_dictionaries("satellites", words[1], words[2])
        elif words[0:2] == ["Simple", "repeats:"]:
            populate_dictionaries("simple", words[2], words[3])
        elif words[0:2] == ["Low", "complexity:"]:
            populate_dictionaries("low_complex", words[2], words[3])
        elif words[0:3] == ["Total", "interspersed", "repeats:"]:
            interspersed_repeats = int (words[3])
    f.close()
    return sequences,total_length,total_length_excl_NX_runs,bases_masked,interspersed_repeats,num_elements,len_occupied



def get_unique_file_set(data_directory):
    unique_file_set = set()
    if data_directory[-1] != r"/":
        data_directory = data_directory + r"/"
    directory_listing = glob.glob(data_directory + r"*.fasta")
    for file in directory_listing:
        candidate_name = file.split("_file")[0]
        unique_file_set.add(candidate_name)
    return unique_file_set



def process_file_set(base_name):
    num_elements = {}
    len_occupied = {}
    for i in indices:
        num_elements[i] = 0
        len_occupied[i] = 0

    sequences = 0
    total_length = 0
    total_length_excl_NX_runs = 0
    bases_masked = 0
    interspersed_repeats = 0


    sample_files = glob.glob(base_name+"*.tbl")
    for subfile in sample_files:
        sequences2,total_length2,total_length_excl_NX_runs2,bases_masked2,interspersed_repeats2,num_elements2,len_occupied2 = process_file(subfile)
        sequences += sequences2
        total_length += total_length2
        total_length_excl_NX_runs += total_length_excl_NX_runs2
        bases_masked += bases_masked2
        interspersed_repeats += interspersed_repeats2
        for i in indices:
            num_elements[i] += num_elements2[i]
            len_occupied[i] += len_occupied2[i]
    return sequences,total_length,total_length_excl_NX_runs,bases_masked,interspersed_repeats,num_elements,len_occupied


def write_reports(output_directory):
    unique_set = get_unique_file_set(output_directory)
    #get the  base_names in unique_set:
    #base_name = "/share/RM_OUTPUT/005_SG_Gale_iCLIP_06_3T3-mT2-3_aFlag_S5Aligned.sortedByCoord.out.fasta"
    for base_name in unique_set:
        print (base_name.split("/")[-1])
        sequences2,total_length2,total_length_excl_NX_runs2,bases_masked2,interspersed_repeats2,num_elements2,len_occupied2 = process_file_set(base_name)
        f = open (base_name.split("/")[-1] + ".summary", "w")
        f.write("============================================================================\n")
        f.write("file name:\t" + base_name.split("/")[-1] + "\n")
        f.write("sequences:\t" + str(sequences2) + "\n")
        f.write("total length:\t" + str(total_length2) + " bp\t(" + str(total_length_excl_NX_runs2) + " bp excl N/X-runs)\n")
        #f.write("bases masked:\t" + str(bases_masked2) + " bp ({0:.2f}%)\n".format(100*bases_masked2/total_length2))
        f.write("bases masked:\t" + str(bases_masked2) + " bp (" + str(round(100*bases_masked2/total_length2,2)) + "%)\n")
        f.write("==============================================================================\n")
        f.write("\t\t\tnumber of elements\tlength\t\t\tpercentage\n")
        f.write("\t\t\telements*\t\toccupied\t\tof sequence\n")
        f.write("------------------------------------------------------------------------------\n")

        f.write("SINEs:\t\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["sines"],len_occupied2["sines"],100*len_occupied2["sines"]/total_length2))
        f.write("\tAlu/B1\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["alus"],len_occupied2["alus"],100*len_occupied2["alus"]/total_length2))
        f.write("\tB2-B4\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["b2b4"],len_occupied2["b2b4"],100*len_occupied2["b2b4"]/total_length2))
        f.write("\tIDs\t\t{0}\t\t\t{1} bp\t\t\t{2:.2f} %\n".format(num_elements2["ids"],len_occupied2["ids"],100*len_occupied2["ids"]/total_length2))
        f.write("\tMIRs\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["mirs"],len_occupied2["mirs"],100*len_occupied2["mirs"]/total_length2))
        f.write("\n")

        f.write("LINEs:\t\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["lines"],len_occupied2["lines"],100*len_occupied2["lines"]/total_length2))
        f.write("\tLINE1\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["line1"],len_occupied2["line1"],100*len_occupied2["line1"]/total_length2))
        f.write("\tLINE2\t\t{0}\t\t\t{1} bp\t\t\t{2:.2f} %\n".format(num_elements2["line2"],len_occupied2["line2"],100*len_occupied2["line2"]/total_length2))
        f.write("\tL3/CR1\t\t{0}\t\t\t{1} bp\t\t\t{2:.2f} %\n".format(num_elements2["l3cr1"],len_occupied2["l3cr1"],100*len_occupied2["l3cr1"]/total_length2))
        f.write("\n")

        f.write("LTR elements:\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["ltr"],len_occupied2["ltr"],100*len_occupied2["ltr"]/total_length2))
        f.write("\tERVL\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["ervl"],len_occupied2["ervl"],100*len_occupied2["ervl"]/total_length2))
        f.write("\tERVL-MaLRs\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["ervl_malrs"],len_occupied2["ervl_malrs"],100*len_occupied2["ervl_malrs"]/total_length2))
        f.write("\tERVL-classI\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["erv_classI"],len_occupied2["erv_classI"],100*len_occupied2["erv_classI"]/total_length2))
        f.write("\tERVL-classII\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["erv_classII"],len_occupied2["erv_classII"],100*len_occupied2["erv_classII"]/total_length2))
        f.write("\n")

        f.write("DNA elements:\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["DNA"],len_occupied2["DNA"],100*len_occupied2["DNA"]/total_length2))
        f.write("\thAT-Charlie\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["hat_charlie"],len_occupied2["hat_charlie"],100*len_occupied2["hat_charlie"]/total_length2))
        f.write("\tTcMar-Tigger\t{0}\t\t\t{1} bp\t\t\t{2:.2f} %\n".format(num_elements2["tcmar"],len_occupied2["tcmar"],100*len_occupied2["tcmar"]/total_length2))
        f.write("\n")

        f.write("Unclassified:\t\t{0}\t\t\t{1} bp\t\t\t{2:.2f} %\n".format(num_elements2["unclassified"],len_occupied2["unclassified"],100*len_occupied2["unclassified"]/total_length2))
        f.write("\n")

        f.write("Total interspersed repeats:\t\t\t{0} bp\t\t{1:.2f} %\n".format(interspersed_repeats2,interspersed_repeats2/total_length2))
        f.write("\n\n")

        f.write("Small RNA:\t\t{0}\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["small_rna"],len_occupied2["small_rna"], 100*len_occupied2["small_rna"]/total_length2))
        f.write("\n")\

        f.write("Satellites:\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["satellites"],len_occupied2["satellites"], 100*len_occupied2["satellites"]/total_length2))
        f.write("Simple repeats:\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["simple"],len_occupied2["simple"], 100*len_occupied2["simple"]/total_length2))
        f.write("Low complexity:\t\t{0}\t\t\t{1} bp\t\t{2:.2f} %\n".format(num_elements2["low_complex"],len_occupied2["low_complex"], 100*len_occupied2["low_complex"]/total_length2))
        f.write("==============================================================================\n")
        f.write("\n\n")

        f.write("\t* most repeats fragmented by insertions or deletions\n")
        f.write("\t  have been counted as one element.\n\n")

        f.write("The query species was assumed to be mus musculus.\n")
        f.write("RepeatMasker Combined Database: Dfam_Consensus-20170127, RepBase-20170127\n\n")
        f.write("Run with rmblastn version 2.6.0+\n")
        f.close()

if len(sys.argv) != 2:
	print ("usage:  rm_report_writer.py  <directory containing tbl files>\n")
else:
        write_reports(sys.argv[1])
