#!/usr/bin/env python

import os
import sys
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file", help="input file; [Default: %(default)s]", type=str, action="store", default="datasets/miniaod/double_muon_low_mass_2018_miniaod.txt")
    parser.add_argument("-t", "--tag", help="tag of the output dataset; [Default: %(default)s]", type=str, action="store", default="")
    parser.add_argument("-j", "--json-file", help="json file; [Default: %(default)s]", type=str, action="store", default="Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt")
    parser.add_argument("-o", "--path-to-store",help="path to user storage; [Default: %(default)s]", type=str, action="store", default="/store/user/bjoshi")
    parser.add_argument("-d", "--data-type", help="data/mc; [Default: %(default)s]", type=str, action="store", default="mc")
    parser.add_argument("-gt", "--global-tag", help="global tag; [Default: %(default)s]", type=str, action="store", default="globalTag=102X_upgrade2018_realistic_v20")
    parser.add_argument("-dt", "--data-tier", help="miniaod/aod; [Default: %(default)s]", type=str, action="store", default="aod")
    args = parser.parse_args()

    templatefile = open("crab_cfg.py",'r')
    text = templatefile.readlines()
    templatefile.close()

    dataset_file = open(args.input_file, 'r')
    datasets = dataset_file.readlines()
    dataset_file.close()
    
    for d in datasets:
       output_tag = args.tag+d.rstrip("\r\n")
       name = output_tag.replace("/","_")
       name = name.replace('-','_')
       filename = 'crab_cfg_'+name[1:]+'.py'
       file_ = open(filename, 'w')
       for line in text:
          if 'requestName' in line: file_.write('config.General.requestName = \''+name+'\'\n')
          elif 'psetName' in line: 
              if (args.data_tier=='miniaod'): file_.write("config.JobType.psetName = \'../muonPogNtuples_miniAOD_cfg.py\'")
              elif (args.data_tier=='aod'): file_.write("config.JobType.psetName = \'../muonPogNtuples_cfg.py\'")
          elif 'workArea' in line: file_.write('config.General.workArea = \''+args.tag+'\'\n')
          elif 'inputDataset' in line: file_.write('config.Data.inputDataset = \''+d.rstrip("\r\n")+'\'\n')
          elif 'outputDatasetTag' in line: file_.write('config.Data.outputDatasetTag = \''+name.replace('-','_')+'\'\n')
          elif 'outLFNDirBase' in line: file_.write('config.Data.outLFNDirBase = \''+args.path_to_store+'\'\n')
          elif 'globalTag' in line: file_.write('config.JobType.pyCfgParams = [\'globalTag='+args.global_tag+'\',\n')
          elif 'runOnMC' in line and args.data_type=='mc': file_.write('                              \'runOnMC=True\',\n')
          elif 'splitting' in line and args.data_type=='mc': file_.write("config.Data.splitting = 'FileBased' \n")
          elif 'unitsPerJob' in line and args.data_type=='mc': file_.write("config.Data.unitsPerJob = 5 \n")
          elif 'lumiMask' in line and args.data_type=='mc': file_.write('\n')
          elif 'runOnMC' in line and args.data_type=='data': file_.write('                              \'runOnMC=False\',\n')
          elif 'splitting' in line and args.data_type=='data': file_.write("config.Data.splitting = 'LumiBased' \n")
          elif 'unitsPerJob' in line and args.data_type=='data': file_.write("config.Data.unitsPerJob = 100 \n")
          elif 'lumiMask' in line and args.data_type=='data': file_.write('config.Data.lumiMask = \''+args.json_file+'\'\n')
          else: file_.write(line)
       
       file_.close()

