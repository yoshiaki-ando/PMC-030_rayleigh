import os
with open("Latitude_Longitude2.dat") as f:
  latlon = f.readlines()

for d in latlon:
  c = d.strip().split(" ")
  print("./main " + c[0] + " " + c[1])
  os.system("./main "+ c[0] + " " + c[1])
