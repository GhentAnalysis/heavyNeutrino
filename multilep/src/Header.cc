#include "heavyNeutrino/multilep/interface/Header.h"

//include c++ library classes
#include <fstream>
#include <sstream>

void Header::addText( const std::string& filePath ){
    std::ifstream textStream( filePath.c_str() );
    std::string line;
    while( std::getline( textStream, line ) ){
        lines.push_back( line );
    }
    textStream.close();
}


Header::Header( const std::vector< std::string >& inputPaths ){
    for( const auto& f : inputPaths ){
        addText( f );
    }
}


void Header::print( std::ostream& os ) const{
    for( const auto& line : lines ){
        os << line << "\n";
    }
    os << std::flush;
}
