#!/usr/bin/env python
# Jeremy Liu
# February 2015
# Adapted from Dan Blackenburg's sample data manager

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

    # Select download URL, file name, data table name, and path using motif_db selector variable
    if motif_db == "pouya":
        BGZ = ['http://compbio.med.harvard.edu/motif-enrichment/pouya_motifs.bed.bgz',
                "pouya_motifs.bed.bgz", "pouya_bgz", "Pouya Encode Motifs (hg19) BGZ"]
        TBI = ['http://compbio.med.harvard.edu/motif-enrichment/pouya_motifs.bed.bgz.tbi',
                "pouya_motifs.bed.bgz.tbi", "pouya_tbi", "Pouya Encode Motifs (hg19) TBI"]
    elif motif_db == "jaspar":
        BGZ = ['http://compbio.med.harvard.edu/motif-enrichment/jaspar_jolma_motifs.bed.bgz',
                "jaspar_jolma_motifs.bed.bgz", "jaspar_bgz", "Jaspar and Jolma Motifs (hg19) BGZ"]
        TBI = ['http://compbio.med.harvard.edu/motif-enrichment/jaspar_jolma_motifs.bed.bgz.tbi',
                "jaspar_jolma_motifs.bed.bgz.tbi", "jaspar_tbi", "Jaspar and Jolma Motifs (hg19) TBI"]
    elif motif_db == "mouse":
        BGZ = ['http://compbio.med.harvard.edu/motif-enrichment/mm9_motifs_split.bed.bgz',
                "mm9_motifs_split.bed.bgz", "mouse_bgz", "Mouse Motifs (mm9) BGZ"]
        TBI = ['http://compbio.med.harvard.edu/motif-enrichment/mm9_motifs_split.bed.bgz.tbi',
                "mm9_motifs_split.bed.bgz.tbi", "mouse_tbi", "Mouse Motifs (mm9) TBI"]
    else:
        BGZ = ['http://compbio.med.harvard.edu/motif-enrichment/pouya_test_motifs.bed.bgz', 
               "pouya_test_motifs.bed.bgz", "test_bgz", "Test Pouya Subset (hg19) BGZ"]
        TBI = ['http://compbio.med.harvard.edu/motif-enrichment/pouya_test_motifs.bed.bgz.tbi',
               "pouya_test_motifs.bed.bgz.tbi", "test_tbi", "Test Pouya Subset (hg19) TBI"]

    # Save and add motif bgz file to motif_databases data table
    bgz_reader = urllib2.urlopen( BGZ[0] )
    bgz_data_table_entry = _stream_fasta_to_file( bgz_reader, target_directory, params,
                            BGZ[1], BGZ[2], BGZ[3] )
    _add_data_table_entry( data_manager_dict, 'motif_databases', bgz_data_table_entry )

    # Save and add motif tbi file to motif_databases data table
    tbi_reader = urllib2.urlopen( TBI[0] )
    tbi_data_table_entry = _stream_fasta_to_file( tbi_reader, target_directory, params,
                            TBI[1], TBI[2], TBI[3] )
    _add_data_table_entry( data_manager_dict, 'motif_databases', tbi_data_table_entry )

def _add_data_table_entry( data_manager_dict, data_table, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][data_table] = data_manager_dict['data_tables'].get( data_table, [] )
    data_manager_dict['data_tables'][data_table].append( data_table_entry )
    return data_manager_dict

def _stream_fasta_to_file( fasta_stream, target_directory, params,
                        fasta_base_filename, value, name, close_stream=True ):
    fasta_filename = os.path.join( target_directory, fasta_base_filename )
    fasta_writer = open( fasta_filename, 'wb+' )
    
    while True:
        buffer = fasta_stream.read(CHUNK_SIZE)
        if not buffer:
            break

        fasta_writer.write(buffer)

    fasta_stream.close()
    fasta_writer.close()
    
    return dict( value=value, name=name, path=fasta_base_filename )

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
    download_motif_databases( data_manager_dict, params, target_directory, options.motif_db )
    
    #save info to json file
    open( filename, 'wb' ).write( to_json_string( data_manager_dict ) )
        
if __name__ == "__main__": main()
