import os 

print(snakemake.input["contamination_summary"])
if len(snakemake.input["contamination_summary"]) > 0:
    file_names = str(snakemake.input["contamination_summary"]).split(' ')
    new_file_names = [cur_file.replace("_contamination_summary.out","") for cur_file in file_names]
    
    f = open(snakemake.params["result_dir"] + "all_contamination_summary.out", "a")
    write_line = ["Sample","Total Reads","Contam Reads", "\n"]
    f.write("\t".join(write_line))                                                             
    
    for i in range(len(file_names)):
        s = open(file_names[i],'r')
        line = s.readlines()
        for cur_line in line:
            file_summary = cur_line.rstrip("\n").split(' ')
            file_summary[1] = file_summary[1] + " (" + str(float(file_summary[1])/float(file_summary[0])) + "%)"
            write_line = [new_file_names[i], file_summary[0], file_summary[1], "\n"]
            f.write("\t".join(write_line))
        s.close()
        
    f.close()
