import glob
import os
import time
import sys
import argparse
import numpy as np

# modify below to match your setup
in_prefix = '/eos/cms/store/group/phys_heavyions/vavladim/'
out_prefix = '/eos/cms/store/group/phys_heavyions/vavladim/condorOutputs/'

parser = argparse.ArgumentParser( description = 'Condor script for generating the .condor file to be passed on to condor to process the input_dataset' )
parser.add_argument( 'input_dataset', help = 'Input dataset to be used. Will go to the directory and loop over the root files.' )
parser.add_argument( 'output_dir', help = 'Output directory.' )
parser.add_argument( 'work_dir', help = 'Location of the working directory.' )
parser.add_argument( 'mode', help = 'Trigger mode.' )

parser.add_argument( '--batches', required = False, default = 20, help = 'Number of batches to divide every 1000 file set.' )

args = parser.parse_args()
in_dataset = args.input_dataset
out_dir = args.output_dir
in_batch = int(args.batches)
workdir = args.work_dir
trig_mode = args.mode

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
script = workdir + '/runtime_wrapper.sh'


out_prefix = out_prefix + in_dataset + '/'

os.system("myschedd bump")
for iS in range(0, len(fp_subsample)):
    if iS == 2:
        break
    prefix_log = workdir + "/logs/" + short_subsample[iS] + "/"
    folder = out_prefix + short_subsample[iS]
    out_dir = folder + '_outdir/'
    log_dir = prefix_log
    filelist = []
    print("Now inside directory ", short_subsample[iS])
    # can have up to 10 subfolders, 1000 jobs each from CRAB
    llist = [[] for x in range(0,10)]
    for path, subdirs, files in os.walk(fp_subsample[iS]):
        for name in files:
            if ".txt" in name:
                print("This file was skipped for being a text file! ", os.path.join(path, name))
            elif "part" in name:
                print("This file was not fully downloaded! ", os.path.join(path, name))
            elif "HiForestMiniAOD" in name:
               # print "Appending filename ",  os.path.join(path, name), " to list"
                for i in range(0,10):
                    if '/000'+str(i)+"/" in os.path.join(path, name):
                        llist[i].append(os.path.join(path, name))
    #create outdir and logdir if they are not already created
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    # index=0
    #open the file which will be passed to condor to put all inputs inside
    # if'/' in short_subsample[iS]:
    extra_file_location = workdir + "/condor_submission/" + in_dataset + '/' + short_subsample[iS] + "/"
    if not os.path.exists(extra_file_location):
        print("Creating directory", extra_file_location)
        os.makedirs(extra_file_location)
    submission_filename = extra_file_location + "submission_script.condor"
    list_files = []
    
    for i in range(0,len(llist)):
        chunked_list = np.array_split(llist[i], in_batch)
        for j in range(0, in_batch):
            textfile = extra_file_location + str(i) + "_subdir_" + str(j) + "_chunk_lists.txt"
            list_files.append(textfile)
            in_list = open(textfile, 'a')
            for item in chunked_list[j]:
                in_list.write(f"\n{item}")
            in_list.close()

    # f = open(submission_filename, 'a')
    # #remove previous contents
    # f.truncate(0)
    # # f.write("MY.SingularityImage     = " + '"/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/batch-team/containers/plusbatch/el8-submit:latest"\n')
    # f.write("MY.WantOS = \"el8\"\n")
    # f.write("Universe   = vanilla\n")
    # f.write("Executable = " + script + "\n")
    # f.write("JobBatchName = \"" + short_subsample[iS] + "\"\n")
    # f.write("+JobFlavour = \"workday\"\n")
    # #idk why we set this to NO and then move the output manually in StepProcess.sh, I just copied it from Leticia's L1 script, it works, if someone feels adventurous, can try to optimise
    # f.write("should_transfer_files = YES\n")
    # f.write("when_to_transfer_output = ON_EXIT\n");
    # f.write("Transfer_Input_Files    = "+ workdir + "/StepProcess.sh," + workdir + "/runtime_wrapper.sh," + workdir + "/runAsymmetry," + workdir + "/libl2residuals.so," + workdir + "/data," + workdir + "/analysis_config.toml")
    # f.write("\n")
    # for in_file in list_files:
    #     only_name = in_file.split('/')[-1]
    #     only_name = only_name.split('.txt')[0]
    #     out_file = out_dir + only_name + '.root'
    #     log_file = log_dir + only_name + '_log.txt'
    #     log_file_e = log_dir + only_name + '_error.txt'
    #     log_file_out = log_dir + only_name + '_out.txt'
    #     job_name = only_name
    #     f.write("Arguments = StepProcess.sh " + workdir + in_file + " " + out_file + " " + trig_mode + "\n")
    #     f.write("Output    = " + log_file_out + "\n")
    #     f.write("Error    = " + log_file_e + "\n")
    #     f.write("Log    = " + log_file + "\n")
    #     f.write("Queue\n")
    # f.close()
    # if iS % 20 == 0:
        # os.system("myschedd bump")
    # os.system("condor_submit " + submission_filename)

All_files = workdir + "/condor_submission/" + in_dataset
All_files_txt = open(All_files + "/_all_files.txt", 'a')
for path, subdirs, files in os.walk(All_files):
    for name in files:
            if ".txt" in name:
                All_files_txt.write(os.path.join(path, name)+"\n")
All_files_txt.close()