// $Id$
#include <string>

#include "stringtools.hpp"

// #include <boost/tokenizer.hpp>

using namespace NRLib2;

std::string
NRLib2::GetPath(const std::string& fName)
{
  static const std::basic_string <char>::size_type npos = -1;

  std::basic_string <char>::size_type index = fName.find_last_of("/");
  if(index == npos)
    index = fName.find_last_of("\\");

  std::string result = "";
  if(index != npos)
    result = fName.substr(0,index+1);

  return(result);
}

std::string
NRLib2::RemovePath(const std::string& fName)
{
  static const std::basic_string <char>::size_type npos = -1;

  std::basic_string <char>::size_type index = fName.find_last_of("/");
  if(index == npos)
    index = fName.find_last_of("\\");

  std::string result = fName;
  if(index != npos)
    result = fName.substr(index+1,fName.length());

  return(result);
}


// std::vector<std::string> NRLib2::GetTokens(const std::string& /*s*/)
//{
//  throw Exception("Not implemented");
//  boost::tokenizer<> tok(s);
//  std::vector<std::string> v;
//  for(boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end();++beg){
//    v.push_back(*beg);
//  }
//  
//  return v;
//}
