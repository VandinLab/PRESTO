import sys
import os

cwd = os.getcwd() + "/"
files_dir = os.path.join(cwd, sys.argv[1])

mapID = {}
ID = 0

for f_ID in range(1,64):
    file_name = "file_" + str(f_ID) 
    file_path = files_dir + file_name
    with open(file_path, 'r') as f:
        for line in f:
            chunks = line.strip().split(" ")
            if len(chunks) != 3:
                print("Something strange")
                continue
            src = chunks[0]
            dst = chunks[1]
            time = chunks[2]
            if src in mapID:
                src_id = mapID[src]
            else:
                src_id = ID
                mapID[src] = src_id
                ID += 1
            if dst in mapID:
                dst_id = mapID[dst]
            else:
                dst_id = ID
                mapID[dst] = dst_id
                ID += 1
            print(str(src_id) + " " + str(dst_id)+ " " + time[5:-3])
#print("Max ID is: ", ID)
