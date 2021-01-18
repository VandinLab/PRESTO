These are the scripts that we used to generate dataset EquinixChicagoA and EquinixChicagoB
the way such code is intended to be used is the following:

REQUIREMENTS:
- Python 3.6 or upper
- Linux-like OS
- Tshark for Python

PROCEDURE:
1. Suppose you douwnloaded the CAIDA dataset, after requesting it at: <link>
2. Locate the folder with the pcap files, note that we separeted dir-a and dir-b in this phase with the 'mv' command, 
   we will thus describe only the procedure for dir-a, for dir-b is the same
3. Place the scripts in the folder containing the folder with the pcap files and execute 'mkdir first_pass' 
4. Launch the script launch_parser.sh with argument the folder containing the 62 pcap files (either dir-a or dir-b)
5. After the procedure is completed (it may take a while), launch 'python dataset_creator.py first_pass/ > dataset_name.txt'
6. Your dataset is in the file 'dataset_name.txt', good work!

NOTE: all the commands are intended to be run in the folder where you placed the scripts, which is specified in point 3.
