/*
 * main.cxx
 *
 *  Created on: Sep 13, 2016
 *      Author: zhangrui
 */
#include "PFlowMonitor.h"

int main(int argc, char* argv[] ) {
  char* option = "files.txt";
  char* folder = "default";
  if( argc > 1 ) {
    option = argv[ 1 ];
  }
  if( argc > 2 ) {
    folder = argv[ 2 ];
  }
  PFlowMonitor monitor;
  std::cout<<option<<std::endl;
  monitor.run(option, folder);
}

