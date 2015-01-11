#!/usr/bin/env python
import optparse

def main():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-d', '--dbkey_description', dest='dbkey_description', action='store', type="string", default=None, help='dbkey_description' )
    (options, args) = parser.parse_args()
    
    filename = args[0]
    print(args)

main()