#include <vector>
#include <string>
#include <Rcpp.h> 
#include "RRA.h"


using namespace Rcpp;
// [[Rcpp::export]]
List rcpp_hello_world() {
    CharacterVector x =
       CharacterVector::create("foo", "bar");
    NumericVector y   =
       NumericVector::create( 0.0, 1.0 ) ;
    //int z0=RRA_main_test();
    List z            = List::create( x, y ) ;
    return z ; 
}
// [[Rcpp::export]]
int call_rra(std::vector<std::string> & RRA_parameters){
  //std::vector<std::string>  RRA_parameters;
  //RRA_parameters.push_back(std::string("RRA"));
  int n=RRA_parameters.size();
  Rcout<<"Calling RRA parameters:\n";
  vector<string> vs;
  for(int i=0; i< RRA_parameters.size();i++){
    Rcout<<RRA_parameters[i]<<" ";
    vs.push_back(RRA_parameters[i]);
  }
  Rcout<<"\n";
  RRA_main(vs);
  return 0;
}
