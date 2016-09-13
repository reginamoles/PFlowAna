/*
 * main.cxx
 *
 *  Created on: Sep 13, 2016
 *      Author: zhangrui
 */
#include "PFlowMonitor.h"

int main(int argc, char* argv[] ) {
  char* option = "files.txt";
  if( argc > 1 ) {
    option = argv[ 1 ];
  }
  PFlowMonitor monitor;
  std::cout<<option<<std::endl;
  monitor.run(option);
}


