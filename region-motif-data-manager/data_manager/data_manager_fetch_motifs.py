#!/usr/bin/env python
#Dan Blankenberg

import sys
import os
import tempfile
import shutil
import optparse
import urllib2
#import uuid
from ftplib import FTP
import tarfile
import zipfile
import gzip
import bz2

from galaxy.util.json import from_json_string, to_json_string

CHUNK_SIZE = 2**20 #1mb

def download_motif_databases( data_manager_dict, params, target_directory, motif_db ):
    TEST_DOWNLOAD_URL = 'http://gehlenborg.com/wp-content/uploads/motif/pouya_test_motifs.bed.bgz'
    #TO DO: Add tbi file

    url = TEST_DOWNLOAD_URL
    fasta_reader = urllib2.urlopen( url )
    
    data_table_entry = _stream_fasta_to_file( fasta_reader, target_directory, params )
    _add_data_table_entry( data_manager_dict, data_table_entry )

def _add_data_table_entry( data_manager_dict, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables']['all_fasta'] = data_manager_dict['data_tables'].get( 'all_fasta', [] )
    data_manager_dict['data_tables']['all_fasta'].append( data_table_entry )
    return data_manager_dict

def _stream_fasta_to_file( fasta_stream, target_directory, params, close_stream=True ):
    fasta_base_filename = "pouya_test_motifs.bed.bgz"
    fasta_filename = os.path.join( target_directory, fasta_base_filename )
    fasta_writer = open( fasta_filename, 'wb+' )
    
    while True:
        buffer = fasta_stream.read(CHUNK_SIZE)
        if not buffer:
            break

        fasta_writer.write(buffer)

    fasta_stream.close()
    fasta_writer.close()
    
    return dict( value="test", name="Test Pouya Subset (hg19)", path=fasta_base_filename )

def main():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-m', '--motif_db', dest='motif_db', action='store', type="string", default=None, help='motif_db' )
    (options, args) = parser.parse_args()
    
    filename = args[0]
    
    params = from_json_string( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}
    
    #Fetch the Motif Database
    download_motif_databases( data_manager_dict, params, target_directory, motif_db )
    
    #save info to json file
    open( filename, 'wb' ).write( to_json_string( data_manager_dict ) )
        
if __name__ == "__main__": main()
