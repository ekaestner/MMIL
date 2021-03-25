#!/usr/bin/env python

import sys, csv

def main(argv):
  if len(sys.argv) != 4:
     print "usage: <control.csv> <column header string> <step>"
     sys.exit(-1)
  f=sys.argv[1]
  c=sys.argv[2]
  s=int(sys.argv[3])

  csvfile = open(f, 'rb')
  reader = csv.reader(csvfile)
  header = reader.next()
  try:
    idx = header.index(c)
  except ValueError:
    print "Error: column does not exist"
    sys.exit(-1)

  all = list( reader )
  try:
    print all[s][idx]
  except IndexError:
    print "Error: step does not exist"
    sys.exit(-1)

if __name__ == "__main__":
  main(sys.argv[1:])