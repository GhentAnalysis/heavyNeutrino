#ifndef Header_H
#define Header_H

//include c++ library classes
#include <string>
#include <vector>
#include <iostream>

class Header{

    public:
        Header( const std::vector< std::string >& );
        
        void print( std::ostream& os = std::cout ) const;

    private:
        void addText( const std::string& filePath );
        std::vector< std::string > lines;

};
#endif
