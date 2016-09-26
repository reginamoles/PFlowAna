/*
 * main.cxx
 *
 *  Created on: Sep 13, 2016
 *      Author: zhangrui
 */
#include "PFlowMonitor.h"
#include "assert.h"

int main(int argc, char* argv[] ) {
  char* option = "files.txt";
  char* folder = "";
  assert(argc > 2);
  option = argv[1];
  folder = argv[2];
  PFlowMonitor monitor;
  std::cout<<option<<std::endl;
  monitor.run(option, folder);
}


