//
// file boost_tuples_and_bind.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This file has a magic spell that I found on the web to make
// boost::bind work with boost::tuples. It defines a function that
// extracts the appropriate element for putting into STL algorithms and stuff

// Example use would be:
// typedef boost::tuples<string,int,int> TUP;
// std::sort( v1.begin() , v1.end() ,
//            bind( less<string>() , bind( &getTupleElement<0,TUP> , _1 ) ,
//                  bind( &getTupleElement<0,TUP> , _2 ) ) );

#include <boost/tuple/tuple.hpp>

template< size_t N, class Tup > inline
typename boost::tuples::element<N,Tup>::type const& getTupleElement( const Tup& t) {

   return boost::tuples::get<N>(t);

}
