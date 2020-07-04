#ifndef RangedMap_H
#define RangedMap_H

/*
Class to map a range of double to another object 
copied from ewkino framework: https://github.com/wverbeke/ewkino/blob/master/Tools/interface/RangedMap.h
*/


//include c++ library classes
#include <map>
#include <initializer_list>
#include <utility>
#include <stdexcept>



template < typename T > class RangedMap {
    
    public:
        using iterator = typename std::map< double, T >::iterator;
        using const_iterator = typename std::map< double, T >::const_iterator;
        //using value_type = typename std::map< double, T >::value_type;
        
        RangedMap( const std::initializer_list< std::pair< const double, T > >& );
        RangedMap( const std::map< double, T >& );
        RangedMap() = default;

        //unlike with an std::map, no new elements can be added with the index operator 
        const T& operator[]( const double ) const;

        iterator begin(){ return lowerBoundMap.begin(); }
        const_iterator begin() const{ return lowerBoundMap.cbegin(); }
        const_iterator cbegin() const{ return lowerBoundMap.cbegin(); }

        iterator end(){ return lowerBoundMap.end(); }
        const_iterator end() const{ return lowerBoundMap.cend(); }
        const_iterator cend() const{ return lowerBoundMap.cend(); }

        bool empty() const{ return lowerBoundMap.empty(); }

    private:
        std::map< double, T > lowerBoundMap;
};



template< typename T > RangedMap< T >::RangedMap( const std::initializer_list< std::pair< const double, T > >& pairList ) :
    lowerBoundMap( pairList )
{}


template< typename T > RangedMap< T >::RangedMap( const std::map< double, T >& boundMap ) :
    lowerBoundMap( boundMap )
{}


template< typename T > const T& RangedMap< T >::operator[]( const double value ) const{

    //make sure the map is not empty
    if( empty() ){
        throw std::out_of_range( "Trying to index empty RangedMap." );
    }

    //throw error if given value is smaller than the first lower bound 
    if( value < (*lowerBoundMap.begin()).first ){
        throw std::invalid_argument( "Given value is smaller than the lower bound of the RangedMap." );
    }

    //look for first map entry not less than given value
    auto it = lowerBoundMap.lower_bound( value );
    
    //the value belongs in the bin below the returned one!    
    --it;
    return (*it).second;
}


#endif
