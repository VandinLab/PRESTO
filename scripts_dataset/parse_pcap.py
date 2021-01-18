import pyshark
import sys
import os

script_dir = os.path.dirname(__file__)
file_path = os.path.join(script_dir, sys.argv[1])
cap = pyshark.FileCapture(file_path, keep_packets=False)

while True: 
        try:
            pack= cap.next()
        except:
            break
        try:
                print(pack['ip'].src + " " + pack['ip'].dst + " " +  pack.sniff_timestamp.replace('.',''))
        except:
                continue
