import glob
import os
import time
import sys
import argparse
import numpy as np

# modify below to match your setup
in_prefix = '/eos/cms/store/group/phys_heavyions/vavladim/'
# in_prefix = '/eos/cms/store/group/phys_heavyions/hbossi/mc_productions/QCD-dijet_pThat15-event-weighted_TuneCP5_5p36TeV_pythia8/'
# in_prefix = '/eos/cms/store/group/phys_heavyions/hbossi/OOJetSubstructure/'

parser = argparse.ArgumentParser( description = 'Condor script for generating the .condor file to be passed on to condor to process the input_dataset' )
parser.add_argument( 'input_dataset', help = 'Input dataset to be used. Will go to the directory and loop over the root files.' )
parser.add_argument( 'text_dir', help = 'Location of list of lists.' )
parser.add_argument( '--batches', required = False, default = 10, help = 'Number of batches to divide every 1000 file set.' )

args = parser.parse_args()
in_dataset = args.input_dataset
in_batch = int(args.batches)
textdir = args.text_dir

in_dir = in_prefix + in_dataset + '/'
subsample_N = len(os.listdir(in_dir))
fp_subsample = [in_dir]
short_subsample = [in_dataset]
# print(fp_subsample)
if subsample_N > 1:
    fp_subsample = os.listdir(in_dir)
    short_subsample = [s for s in fp_subsample]
    fp_subsample = [in_dir + s for s in fp_subsample]
print("Submitting subsample:", short_subsample)

# os.system("myschedd bump")
for iS in range(0, len(fp_subsample)):
    # if iS == 2:
    #     break
    filelist = []
    print("Now inside directory ", short_subsample[iS])
    # print("Now inside directory ", fp_subsample[iS])
    # can have up to 10 subfolders, 1000 jobs each from CRAB
    llist = [[] for x in range(0,10)]
    for path, subdirs, files in os.walk(fp_subsample[iS]):
        for name in files:
            # print("Name that gets checked", name)
            if ".txt" in name:
                print("This file was skipped for being a text file! ", os.path.join(path, name))
            elif "part" in name:
                print("This file was not fully downloaded! ", os.path.join(path, name))
            elif ".root" in name:
               # print "Appending filename ",  os.path.join(path, name), " to list"
                for i in range(0,10):
                    if '/000'+str(i)+"/" in os.path.join(path, name):
                        llist[i].append(os.path.join(path, name))
                        # print(llist[i])
    #create outdir and logdir if they are not already created
    # index=0
    #open the file which will be passed to condor to put all inputs inside
    # if'/' in short_subsample[iS]:
    extra_file_location = textdir + "/" + in_dataset + "/" + short_subsample[iS] + "/"
    if not os.path.exists(extra_file_location):
        print("Creating directory", extra_file_location)
        os.makedirs(extra_file_location)
    list_files = []
    
    for i in range(0,len(llist)):
        chunked_list = np.array_split(llist[i], in_batch)
        for j in range(0, in_batch):
            textfile = extra_file_location + str(i) + "_subdir_" + str(j) + "_chunk_lists.txt"
            list_files.append(textfile)
            in_list = open(textfile, 'a')
            in_list.truncate(0)
            for item in chunked_list[j]:
                in_list.write(f"\n{item}")
            in_list.close()

All_files = textdir + "/" + in_dataset + "/"
All_files_txt = open("data/txt/" + in_dataset + ".txt", 'a')
print("All files contained in list ", "data/txt/" + in_dataset + ".txt")
All_files_txt.truncate(0)
for path, subdirs, files in os.walk(All_files):
    for name in files:
            if ".txt" in name:
                All_files_txt.write(os.path.join(path, name)+"\n")
All_files_txt.close()