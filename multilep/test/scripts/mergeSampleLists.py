import sys
import os


def buildSampleList( input_file_name ):
    sampleList = []
    with open( input_file_name ) as f:
        for line in f.readlines():
            line = line.strip()

            #remove empty lines or comments
            if line.startswith('#') or not line:
                continue

            #remove option lines
            if line.startswith('+') or line.startswith('-'):
                continue

            #remove skim
            line = line.split(':')[-1]

            sampleList.append( line )

    #check for duplicate samples (i.e faulty sample lists)
    checkSampleListConsistency( sampleList )

    return sampleList


def checkSampleListConsistency( sample_list ):
    sample_set = set( sample_list )
    if len( sample_set ) != len( sample_list ):
        raise AssertionError( 'Duplicate samples are present!' )


#sorting key for sample list sorting
def sampleSortKey( sample_name ):
    
    #data will be put on top in the merged sample list
    if not 'SIM' in sample_name :
        return ( 0, sample_name )
    else:
        return ( 1, sample_name )


def mergeSampleLists( lhs, rhs ):
    merged_list = lhs + rhs 

    #remove overlapping samples
    merged_set = set( merged_list )
    merged_list = list( merged_set )

    #sort the samples to make the list cleaner
    merged_list = sorted( merged_list, key = sampleSortKey )

    return merged_list


def writeSampleList( sample_list, skim, output_file_name ):
    if os.path.isfile( output_file_name ):
        raise IOError( "File '{}' already exists! Overwiting sample list is not allowed.".format( output_file_name ) )

    with open( output_file_name, 'w' ) as f:
        for sample in sample_list:
            f.write( '{}:{}\n'.format( skim, sample ) )
            


if __name__ == '__main__':
    if len( sys.argv ) != 5 :
        error_message = 'This script requires 4 command line arguments to run.\n'
        error_message += 'Usage : python mergeSampleLists.py <sample_list_1> <sample_list_2> <skim> <new_sample_list_name>'
        raise RuntimeError( error_message ) 
    else :
        lhs_sample_list_file_name = sys.argv[1] 
        if not os.path.isfile( lhs_sample_list_file_name ):
            raise IOError( "Specified sample list '{}' does not exist.".format( lhs_sample_list_file_name ) )
        rhs_sample_list_file_name = sys.argv[2]
        if not os.path.isfile( rhs_sample_list_file_name ):
            raise IOError( "Specified sample list '{}' does not exist.".format( rhs_sample_list_file_name ) )

        lhs_sample_list = buildSampleList( lhs_sample_list_file_name )
        rhs_sample_list = buildSampleList( rhs_sample_list_file_name )
        merged_sample_list = mergeSampleLists( lhs_sample_list, rhs_sample_list )

        skim = sys.argv[3] 
        new_sample_list_file_name = sys.argv[4]
        writeSampleList( merged_sample_list, skim, new_sample_list_file_name )
